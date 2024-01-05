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

include {SAMTOOLS_SAMPLES_01} from '../samtools/samtools.samples.01.nf'
include {moduleLoad;isHg19;isHg38;getVersionCmd;getGnomadGenomePath} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_CONCAT_01} from '../bcftools/bcftools.concat.01.nf'
include {JVARKIT_VCF_TO_INTERVALS_01} from '../jvarkit/jvarkit.vcf2intervals.nf'

workflow ULTRA_RARES_01 {
	take:
		meta
		reference
		vcf
		bams
		bed
	main:
		version_ch = Channel.empty()

		bams_ch = SAMTOOLS_SAMPLES_01(["with_header":false,"allow_multiple_references":false,"allow_duplicate_samples":false],reference,file("NO_FILE"),bams)
		version_ch = version_ch.mix(bams_ch.version)		

		ch2 = CHECK_NO_COMMON_SAMPLES_AND_SPLIT_BAMS(meta,vcf,bams_ch.output)
		version_ch = version_ch.mix(ch2.version)
		

		if(bed.name.equals("NO_FILE")) {
			ch3 = JVARKIT_VCF_TO_INTERVALS_01(meta, reference, vcf, file("NO_FILE"))
			version_ch = version_ch.mix(ch3.version)
			intervals = ch3.bed
		} else {
			intervals = Channel.fromPath(bed)
			}
		intervals = intervals.splitCsv(header:false,sep:'\t').
			map{T->T[0]+":"+((T[1] as int)+1)+"-"+T[2]}

		rares_ch = APPLY_ULTRARES(meta, reference, vcf, intervals, ch2.output)
		version_ch = version_ch.mix(rares_ch.version)

		ch5_ch = COLLECT_TO_FILE_01(meta, rares_ch.vcf.collect())
		version_ch = version_ch.mix(ch5_ch.version)

		ch6_ch = BCFTOOLS_CONCAT_01(meta,ch5_ch.output)
		version_ch = version_ch.mix(ch6_ch.version)

		version_ch = MERGE_VERSION(meta, "UltraRares", "utra rares", version_ch.collect())
	emit:
		version = version_ch
		vcf = ch6_ch.vcf
	}

process CHECK_NO_COMMON_SAMPLES_AND_SPLIT_BAMS {
executor "local"
input:
	val(meta)
	val(vcf)
	val(bams)
output:
	path("bams.list"),emit:output
	path("version.xml"),emit:version
script:
	def nbams = meta.n_bams_per_hc_call?:10
"""
hostname 1>&2
${moduleLoad("bcftools")}

bcftools query -l "${vcf}" | sort | uniq > jeter.a
test -s jeter.a

cut -f 1 "${bams}" | sort | uniq > jeter.b
test -s jeter.b

comm -12 jeter.a jeter.b > jeter.c
test ! -s jeter.c

rm jeter.a jeter.b jeter.c

mkdir BAMS
# increasing number of bams per list until we reach nbams
cut -f 3 "${bams}" |  awk 'BEGIN {MAX=${nbams};N=0;M=1;IDX=1;} {N++;OUT=sprintf("BAMS/cluster.%03d.list",IDX);print >> OUT;if(N==M){close(OUT);N=0;IDX++;M++;if(M>MAX) M=MAX;}} END{ close(OUT); }'
find \${PWD}/BAMS -type f -name "cluster*.list" | sort -V > bams.list
test -s bams.list

#####
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">check there is no common sample between bam and vcf</entry>
        <entry key="vcf">${vcf}</entry>
        <entry key="bams">${bams}</entry>
        <entry key="versions">${getVersionCmd("bcftools")}</entry>
</properties>
EOF
"""
}

process APPLY_ULTRARES {
tag "${interval}"
afterScript "rm -rf TMP"
memory "5g"
cpus 4
input:
	val(meta)
	val(reference)
	val(vcf)
	val(interval)
	val(bams)
output:
	path("genotyped.bcf"),emit:vcf
	path("version.xml"),emit:version
script:
	def ploidy = 2
	def dbsnp = ""
	def gnomadGenome= getGnomadGenomePath(meta,reference);
	def gnomadPop = meta.gnomad_population?:"AF_nfe"
	def gnomadAF = meta.gnomad_max_af?:0.01
	def extraBcftoolsView1 = meta.extraBcftoolsView1?:""
	def extraBcftoolsView2 = meta.extraBcftoolsView2?:""
"""
hostname 1>&2
${moduleLoad("bcftools bedtools gatk4 jvarkit")}
set -x
mkdir -p TMP

bcftools query -l "${vcf}" > TMP/samples.txt

##### prepare code filter 

echo -n "final Set<String> CASES = new HashSet<>(Arrays.asList(" > TMP/att.code
sed 's/^/"/;s/\$/"/' TMP/samples.txt | paste -s -d ',' >> TMP/att.code
echo "));" >> TMP/att.code

echo -n "final int n_controls = " >> TMP/att.code
xargs -a "${bams}" cat | wc -l  >> TMP/att.code
echo ";" >> TMP/att.code

cat << EOF >> TMP/att.code

/* overlapping variants */
if(!variant.hasAttribute("ORIGINAL")) return false;

final Set<String> others = new HashSet<>(variant.getAttributeAsStringList("OSAMPLES",""));

others.addAll(variant.getGenotypes().stream().
		filter(G->!(G.isNoCall() ||  G.isHomRef())).
		map(G->G.getSampleName()).
		filter(S->!CASES.contains(S)).
		collect(Collectors.toSet()));

if (others.isEmpty() ) return true;

if ( (others.size()/(double)n_controls) > ${gnomadAF} ) return false;
return new VariantContextBuilder(variant).
	attribute("OSAMPLES",new ArrayList<String>(others)).
	make();
EOF

bcftools view ${extraBcftoolsView1} --regions "${interval}" -O u "${vcf}" |\
	bcftools norm --no-version --multiallelics -any --fasta-ref "${reference}" -O v |\
	awk -F '\t' '/^##/ {print;next;} /^#CHROM/ {printf("##INFO=<ID=OSAMPLES,Number=.,Type=String,Description=\\\"Found in other samples\\\">\\n");printf("##INFO=<ID=ORIGINAL,Number=0,Type=Flag,Description=\\\"PRESENT IN ORIGNAL VCF\\\">\\n");print;next;} {if(\$5=="*") next;OFS="\t";\$8=sprintf("ORIGINAL;%s",\$8);print;}' |\
	bcftools sort --max-mem ${task.memory.giga}G  -T TMP -O z -o TMP/jeter.vcf.gz

## gnomad
if ${!gnomadGenome.isEmpty() && !gnomadPop.isEmpty()} ; then
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/vcfgnomad.jar \
			--bufferSize 10000 \
			--max-af ${gnomadAF} \
			--gnomad "${gnomadGenome}" --fields "${gnomadPop}" TMP/jeter.vcf.gz  |\
	bcftools view -e 'FILTER ~ "GNOMAD_GENOME_BAD_AF"' -O z -o  TMP/jeter2.vcf.gz
	mv TMP/jeter2.vcf.gz TMP/jeter.vcf.gz
fi

if ${!extraBcftoolsView2.isEmpty()} ; then
	bcftools view ${extraBcftoolsView2} -O z -o TMP/jeter2.vcf.gz TMP/jeter.vcf.gz
	mv TMP/jeter2.vcf.gz TMP/jeter.vcf.gz
fi

bcftools index --tbi TMP/jeter.vcf.gz

cat "${bams}" | while read BLIST
do
	SECONDS=0
	echo "bams:\${BLIST}" 1>&2
	echo -n "\${BLIST} " >> TMP/processed.txt
	bcftools query -f '.' TMP/jeter.vcf.gz | wc -c  >> TMP/processed.txt

	# break the loop if no more variant
	if [[ \$(bcftools query -f '.' TMP/jeter.vcf.gz | wc -c ) -le 0 ]] ; then
		break
	fi

	gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" HaplotypeCaller \
		-R "${reference}" \
		${dbsnp.isEmpty()?"":"--dbsnp \"${dbsnp}\" "} \
		${meta.mapq && ((meta.mapq as int) >= 0)?"--minimum-mapping-quality ${meta.mapq}":""} \
		--sample-ploidy "${ploidy}" \
		--do-not-run-physical-phasing \
		--alleles TMP/jeter.vcf.gz \
		-L TMP/jeter.vcf.gz \
		-I "\${BLIST}" \
		-O TMP/jeter2.vcf.gz
	
	bcftools norm -O v --no-version --multiallelics -any --regions "${interval}" --fasta-ref "${reference}" TMP/jeter2.vcf.gz |\
		awk -F '\t' '/^##/ {print;next;} /^#CHROM/ {printf("##INFO=<ID=OSAMPLES,Number=.,Type=String,Description=\\\"Found in other samples\\\">\\n");printf("##INFO=<ID=ORIGINAL,Number=0,Type=Flag,Description=\\\"PRESENT IN ORIGNAL VCF\\\">\\n");print;next;} {if(\$5=="*") next;print;}' |\
		bcftools sort --max-mem ${task.memory.giga}G  -T TMP -O z -o TMP/jeter3.vcf.gz
	mv TMP/jeter3.vcf.gz TMP/jeter2.vcf.gz
	bcftools index -t -f TMP/jeter2.vcf.gz


	bcftools merge -O v TMP/jeter.vcf.gz TMP/jeter2.vcf.gz |\
	java  -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/vcffilterjdk.jar --nocode -f TMP/att.code |\
	bcftools view -O u --samples-file TMP/samples.txt |\
	bcftools view -i 'AC>0' -O z -o TMP/jeter3.vcf.gz

	mv TMP/jeter3.vcf.gz TMP/jeter.vcf.gz
	bcftools index -t -f TMP/jeter.vcf.gz
	
	echo "That took: \${SECONDS} seconds." 1>&2
done

bcftools sort --max-mem ${task.memory.giga}G  -T . -O b -o genotyped.bcf TMP/jeter.vcf.gz
bcftools index genotyped.bcf


####
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">extract rares variants</entry>
        <entry key="vcf">${vcf}</entry>
        <entry key="interval">${interval}</entry>
        <entry key="bams">${bams}</entry>
        <entry key="gnomad.genome">${gnomadGenome}</entry>
        <entry key="gnomad.pop">${gnomadPop}</entry>
        <entry key="gnomad.af">${gnomadAF}</entry>
        <entry key="versions">${getVersionCmd("bcftools bedtools gatk jvarkit/vcffilterjdk jvarkit/vcfgnomad")}</entry>
</properties>
EOF
"""
}
