/*

Copyright (c) 2023 Pierre Lindenbaum

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

def gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults().putReference()


gazoduc.build("vcf", "NO_FILE").
	desc("Indexed VCF file.").
	existingFile().
	required().
	put()

gazoduc.make("pedigree","NO_FILE").
	description("pedigree containing trios").
	existingFile().
	required().
	put()

gazoduc.make("gtf","NO_FILE").
	description("gtf file indexed with tabix").
	existingFile().
	required().
	put()


gazoduc.make("gnomadAF",0.01).
	description("Frequency in gnomad").
	setDouble().
	put()

gazoduc.make("gnomadPop","AF_NFE").
	description("Population in gnomad").
	put()


include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {getVersionCmd;getGnomadGenomePath;runOnComplete;moduleLoad;isBlank;isHg38;isHg19} from '../../modules/utils/functions.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {SIMPLE_PUBLISH_01} from '../../modules/utils/publish.simple.01.nf'


if( params.help ) {
    gazoduc.usage().
	name("HET Composite").
	desc("Scan for het composites").
	print();
    exit 0
} else {
   gazoduc.validate();
}



workflow {
	ch = SCAN_HET_COMPOSITES(params, params.reference, file(params.vcf), file(params.pedigree), file(params.gtf) )
	html = VERSION_TO_HTML(params,ch.version)
	SIMPLE_PUBLISH_01(params, Channel.empty().mix(html.html).mix(ch.version).mix(ch.zip).collect())
	}

runOnComplete(workflow);

workflow SCAN_HET_COMPOSITES {
    take:
	    meta
	    reference
	    vcf
	    pedigree
	    gtf
    main:
		version_ch = Channel.empty()

		ch1 = VCF_TO_BED(meta, vcf)
		version_ch = version_ch.mix(ch1.version)

		ch4= INTER_VCF_GTF(meta, ch1.bed, gtf)
		version_ch = version_ch.mix(ch4.version)

		ctg_vcf_gtf = ch4.output.splitCsv(header:false,sep:'\t')

		ch2 = PER_CONTIG(meta, reference, pedigree, ctg_vcf_gtf)
		version_ch = version_ch.mix(ch2.version)

		ch3 = MERGE(meta, ch2.vcf.collect(), ch2.genes_report.collect(), ch2.variants_report.collect())
		version_ch = version_ch.mix(ch3.version)

		version_ch = MERGE_VERSION(meta, "Het Composites", "Het Composites", version_ch.collect())
    emit:
	    version = version_ch
	    zip = ch3.zip
    }


process INTER_VCF_GTF {
executor "local"
input:
	val(meta)
	path(bed)
	val(gtf)
output:
	path("join.tsv"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("htslib")}

join -t '\t' -1 1 -2 1 -o '2.1,2.2' \
	<( tabix -l "${gtf}" | sort) \
	<( cut -f1,4 "${bed}" | sort -t '\t' -k1,1 -k2,2) |\
	awk '{printf("%s\t${gtf}\\n",\$0);}' > join.tsv

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
	val(reference)
	val(pedigree)
	tuple val(contig), val(vcf), val(gtf)
output:
	path("contig.bcf"),emit:vcf
	path("genes.report"),emit:genes_report
	path("variants.report"),emit:variants_report
	path("version.xml"),emit:version
script:
       	def snpeffdb = (isHg19(reference)?"GRCh37.75":(isHg38(reference)?"GRCh38.86":"TODO"))
        def gnomadVcf = getGnomadGenomePath(meta,reference)
        def gnomadAF = meta.gnomadAF?:0.01
        def gnomadPop = meta.gnomadPop?:"AF_nfe"
        def soacn = meta.soacn?:"SO:0001629,SO:0001818"
"""
hostname 1>&2
${moduleLoad("bcftools jvarkit bedtools snpEff")}

mkdir -p TMP
set -x

awk -F '\t' 'BEGIN {printf("Genotype f,m,c;");} {P=\$6; if((P=="case" || P=="affected") && \$3!="0" && \$4!="0") {printf("c=variant.getGenotype(\\"%s\\");f=variant.getGenotype(\\"%s\\");m=variant.getGenotype(\\"%s\\");if(!c.isHet()) return false; if(f.isHet() && m.isHet()) return false;if(!(f.isHet() || m.isHet())) return false;",\$2,\$3,\$4);} } END {printf("return variant.getGenotypes().stream().noneMatch(G->G.isHomVar());\\n");}'  '${pedigree}' > TMP/jeter.code

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
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffilterjdk -f TMP/jeter.code > TMP/jeter2.vcf
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

bcftools view  -e 'FILTER ~ "GNOMAD_GENOME_BAD_AF" || FILTER ~ "GNOMAD_GENOME_AS_VQSR"' -O v -o TMP/jeter2.vcf TMP/jeter1.vcf
mv TMP/jeter2.vcf TMP/jeter1.vcf

bcftools query -f . TMP/jeter1.vcf |wc -c 1>&2


java -Xmx${task.memory.giga}G  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfcomposite \
	--extractors ANN/GeneId \
	--filter "" \
	--genes genes.report \
	--pedigree  "${pedigree}" \
	--report variants.report \
	--tmpDir TMP \
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
	path("${meta.prefix?:""}archive.zip"),emit:zip
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}

mkdir -p "TMP"

bcftools concat -O b -o "TMP/${meta.prefix?:""}concat.bcf" ${L1.join(" ")}
bcftools index "TMP/${meta.prefix?:""}concat.bcf"

cat ${L2.join(" ")} | LC_ALL=C sort -T . | uniq > "TMP/${meta.prefix?:""}genes.txt"
cat ${L3.join(" ")} | LC_ALL=C sort -T . | uniq > "TMP/${meta.prefix?:""}variants.txt"

mv TMP "${meta.prefix?:""}archive"

zip -r -9  "${meta.prefix?:""}archive.zip"  "${meta.prefix?:""}archive"

cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">merge data</entry>
</properties>
EOF
"""
}
