/*

Copyright (c) 2024 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
The MIT License (MIT)
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
nextflow.enable.dsl=2

include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {dumpParams;getVersionCmd;runOnComplete;moduleLoad;isBlank} from '../../modules/utils/functions.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'

if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}


workflow {
	ch = SCAN_HET_COMPOSITES([:], params.genomeId, file(params.vcf), file(params.pedigree) )
	html = VERSION_TO_HTML(ch.version)
	}

runOnComplete(workflow);

workflow SCAN_HET_COMPOSITES {
    take:
	    meta
	    genomeId
	    vcf
	    pedigree
    main:
		version_ch = Channel.empty()

		ch1 = VCF_TO_BED([with_header:false], vcf)
		version_ch = version_ch.mix(ch1.version)

		ch4= INTER_VCF_GTF([:], genomeId, ch1.bed)
		version_ch = version_ch.mix(ch4.version)

		ctg_vcf_gtf = ch4.output.splitCsv(header:false,sep:'\t')

		ch2 = PER_CONTIG([:], genomeId, pedigree, ctg_vcf_gtf)
		version_ch = version_ch.mix(ch2.version)

		ch3 = MERGE([:], ch2.vcf.collect(), ch2.genes_report.collect(), ch2.variants_report.collect())
		version_ch = version_ch.mix(ch3.version)

		version_ch = MERGE_VERSION("Het Composites", version_ch.collect())
    emit:
	    version = version_ch
	    zip = ch3.zip
    }


process INTER_VCF_GTF {
executor "local"
input:
	val(meta)
	val(genomeId)
	path(bed)
output:
	path("join.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def gtf = params.genomes[genomeId].gtf
"""
hostname 1>&2
${moduleLoad("htslib")}

mkdir -p TMP

join -t '\t' -1 1 -2 1 -o '2.1,2.2' \
	<( tabix -l "${gtf}" | sort -T TMP ) \
	<( cut -f1,4 "${bed}" | sort -T TMP -t '\t' -k1,1 -k2,2) |\
	awk -F '\t' '{printf("%s\t${gtf}\\n",\$0);}' > TMP/join.tsv

mv TMP/join.tsv ./

test -s join.tsv

cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">join vcf and gtf on contigs</entry>
</properties>
EOF
"""
}


process PER_CONTIG {
tag "${contig} ${vcf}"
memory "3g"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(genomeId)
	val(pedigree)
	tuple val(contig), val(vcf), val(gtf)
output:
	path("contig.bcf"),emit:vcf
	path("genes.report"),emit:genes_report
	path("variants.report"),emit:variants_report
	path("version.xml"),emit:version
script:
	def genome =  params.genomes[genomeId]
	def reference = genome.fasta
       	def snpeffdb = genome.snpeff_database_name
        def gnomadVcf = genome.gnomad_genome
        def gnomadAF = params.gnomadAF
        def gnomadPop = params.gnomadPop
        def soacn = params.soacn?:"SO:0001629,SO:0001818"
"""
hostname 1>&2
${moduleLoad("bcftools bedtools snpEff jvarkit")}

mkdir -p TMP
set -x

which java 1>&2
echo "\${JAVA_HOME}" 1>&2

cat << EOF > TMP/jeter.code


private boolean acceptControl(final VariantContext vc,String sm) {
	final Genotype g = vc.getGenotype(sm);
	if(g!=null && g.isHomVar()) return false;
	return true;
	}



private boolean acceptTrio(final VariantContext vc,String cm,String fm,String mm) {
	final Genotype c = vc.getGenotype(cm);
        if(c==null) return false;
	final Genotype m = vc.getGenotype(mm);
	final Genotype f = vc.getGenotype(fm);
	// child must be HET
	if(!c.isHet()) return false;
	// discard de novo
	if(f.isHomRef() && m.isHomRef()) return false;
	// parent shouldn't be homvar
	if(f.isHomVar()) return false;
	if(m.isHomVar()) return false;
	// at least one parent should be het
	if(!(f.isHet() || m.isHet())) return false;
	return true;
	}

public Object apply(final VariantContext variant) {
EOF

## all other samples are controls
comm -13 \
	<(cut -f 2 '${pedigree}' | sort | uniq) \
	<(bcftools query -l '${vcf}'| sort | uniq) |\
	awk '{printf("if(!acceptControl(variant,\\"%s\\")) return false;\\n",\$1);}' >> TMP/jeter.code

awk -F '\t' '(\$6=="control" || \$6=="unaffected") {printf("if(!acceptControl(variant,\\"%s\\")) return false;\\n",\$2);}'  '${pedigree}' >> TMP/jeter.code

awk -F '\t' '((\$6=="case" || \$6=="affected") && \$3!="0" && \$4!="0") {printf("if(acceptTrio(variant,\\"%s\\",\\"%s\\",\\"%s\\")) return true;\\n",\$2,\$3,\$4);}'  '${pedigree}'  >> TMP/jeter.code


echo "return false;}" >> TMP/jeter.code


tabix "${gtf}" "${contig}" |\
	grep -wF protein_coding |\
	awk -F '\t' '(\$3=="exon") {printf("%s\t%d\t%s\\n",\$1,int(\$4)-1,\$5);}' |\
	bedtools slop -i - -g "${reference}.fai" -b 5 |\
	sort -T TMP -t '\t' -k1,1 -k2,2n |\
	bedtools merge > TMP/exons.bed


if ! test -s TMP/exons.bed ; then
	echo -e '${contig}\t0\t1' > TMP/exons.bed
fi


bcftools view --regions-file TMP/exons.bed "${vcf}"  |\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffilterjdk --body -f TMP/jeter.code > TMP/jeter2.vcf
mv TMP/jeter2.vcf TMP/jeter1.vcf

bcftools query -f . TMP/jeter1.vcf |wc -c 1>&2


java -jar -Xmx${task.memory.giga}G  -Djava.io.tmpdir=TMP \${SNPEFF_JAR} eff \
        -config \${SNPEFF_CONFIG} \
        -nodownload -noNextProt -noMotif -noInteraction -noLog -noStats -chr chr -i vcf -o vcf "${snpeffdb}" TMP/jeter1.vcf > TMP/jeter2.vcf
mv TMP/jeter2.vcf TMP/jeter1.vcf

bcftools query -f . TMP/jeter1.vcf |wc -c 1>&2


java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffilterso \
                --acn "${soacn}" TMP/jeter1.vcf > TMP/jeter2.vcf
mv TMP/jeter2.vcf TMP/jeter1.vcf

bcftools query -f . TMP/jeter1.vcf |wc -c 1>&2


java -Xmx${task.memory.giga}G  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfgnomad  --bufferSize 10000 \
                --gnomad "${gnomadVcf}" \
                --fields "${gnomadPop}" \
                --max-af "${gnomadAF}" TMP/jeter1.vcf > TMP/jeter2.vcf
mv TMP/jeter2.vcf TMP/jeter1.vcf

## bcftools ne marche pas avec   -e 'FILTER ~ "GNOMAD_GENOME_BAD_AF" || FILTER ~ "GNOMAD_GENOME_AS_VQSR"' filtre pas forcement 
awk  '\$0 ~ /^#/ || !(\$7 ~ /GNOMAD_GENOME_BAD/)' TMP/jeter1.vcf > TMP/jeter2.vcf
mv TMP/jeter2.vcf TMP/jeter1.vcf

bcftools query -f . TMP/jeter1.vcf |wc -c 1>&2


java -Xmx${task.memory.giga}G  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfcomposite \
	--extractors ANN/GeneId \
	--filter "" \
	--genes genes.report \
	--pedigree  "${pedigree}" \
	--report variants.report \
	--tmpDir TMP \
	--max-variants 30 \
	TMP/jeter1.vcf > TMP/jeter2.vcf

mv TMP/jeter2.vcf TMP/jeter1.vcf


bcftools sort -T TMP -o contig.bcf -O b TMP/jeter1.vcf
bcftools index contig.bcf

cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">get het composite for contig</entry>
        <entry key="versions">${getVersionCmd("bcftools bedtools")}</entry>
</properties>
EOF
"""
}


process MERGE {
tag "N=${L1.size()}"
input:
	val(meta)
	val(L1)
	val(L2)
	val(L3)
output:
	path("${params.prefix?:""}archive.zip"),emit:zip
	path("version.xml"),emit:version
script:
	def prefix = params.prefix?:""
"""
hostname 1>&2
${moduleLoad("bcftools")}

mkdir -p "TMP"

bcftools concat -a -O z -o "TMP/${prefix}concat.vcf.gz" ${L1.join(" ")}
bcftools index "TMP/${prefix}concat.vcf.gz"

cat ${L2.join(" ")} | LC_ALL=C sort -T . | uniq > "TMP/${prefix}genes.txt"
cat ${L3.join(" ")} | LC_ALL=C sort -T . | uniq > "TMP/${prefix}variants.txt"

mv TMP "${prefix?:""}archive"

zip -r -9  "${prefix}archive.zip"  "${prefix}archive"

cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">merge data</entry>
</properties>
EOF
"""
}
