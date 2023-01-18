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
include {getVersionCmd;moduleLoad;isHg19;getGnomadGenomePath} from '../../modules/utils/functions.nf'
include {JVARKIT_VCF_TO_BED_01} from '../jvarkit/jvarkit.vcf2bed.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_CONCAT_01} from '../bcftools/bcftools.concat.01.nf'


workflow FAST_ANNOT_01 {
	take:
		meta
		reference
		vcf
		bed
		samples
	main:
		version_ch = Channel.empty()

		ch3 = JVARKIT_VCF_TO_BED_01(meta, reference, vcf, bed)
		version_ch = version_ch.mix(ch3.version)


		annot_ch = ANNOT(meta, reference, samples, ch3.output.splitCsv(header:false,sep:'\t'))
		version_ch = version_ch.mix(annot_ch.version)

		ch5_ch = COLLECT_TO_FILE_01(meta, annot_ch.vcf.collect())
		version_ch = version_ch.mix(ch5_ch.version)


		ch6_ch = BCFTOOLS_CONCAT_01(meta,ch5_ch.output)
		version_ch = version_ch.mix(ch6_ch.version)

		version_ch = MERGE_VERSION(meta, "FastAnnotVcf", "fast annot vcf", version_ch.collect())
	emit:
		version = version_ch
		vcf = ch6_ch.vcf
	}



process ANNOT {
tag "${bed}"
afterScript "rm -rf TMP"
memory "5g"
input:
	val(meta)
	val(reference)
	path(samples)
	tuple val(bed),val(vcf)
output:
	path("annot.bcf"),emit:vcf
	path("version.xml"),emit:version
script:
	def snpeffdb = (isHg19(reference)?"GRCh37.75":"TODO")
	def gnomadVcf = getGnomadGenomePath(meta,reference)
	def gnomadAF = meta.gnomadAF?:0.01
	def gnomadPop = meta.gnomadPop?:"AF_nfe"
	def soacn = meta.soacn?:"SO:0001629,SO:0001818"
"""
hostname 1>&2
${moduleLoad("bcftools jvarkit snpEff/0.0.0")}
mkdir -p TMP

if ${!samples.name.equals("NO_FILE")} ; then
	bcftools query -l "${vcf}" | sort -T TMP | uniq  > TMP/jeter.a
	cat "${samples}" | sort -T TMP | uniq > TMP/jeter.b
	comm -12 TMP/jeter.a TMP/jeter.b > TMP/jeter.samples
	test -s TMP/jeter.samples
fi

bcftools view --regions-file "${bed}" ${samples.name.equals("NO_FILE")?"":"--samples-file TMP/jeter.samples"} -O b -o TMP/jeter1.bcf "${vcf}"

if ${!samples.name.equals("NO_FILE")} ; then
	bcftools view --min-ac 1 -O b -o TMP/jeter2.bcf TMP/jeter1.bcf
	mv TMP/jeter2.bcf TMP/jeter1.bcf
fi

bcftools view TMP/jeter1.bcf |\
java -jar -Xmx${task.memory.giga}G  -Djava.io.tmpdir=TMP \${SNPEFF_JAR} eff \
	-config \${SNPEFF_CONFIG} \
	-nodownload -noNextProt -noMotif -noInteraction -noLog -noStats -chr chr -i vcf -o vcf "${snpeffdb}" > TMP/jeter1.vcf


if ${!soacn.isEmpty()} ; then
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/vcffilterso.jar \
		--acn "${soacn}" TMP/jeter1.vcf > TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf
fi

java -Xmx${task.memory.giga}G  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/vcfgnomad.jar  --bufferSize 10000 \
		--gnomad "${gnomadVcf}" \
		--fields "${gnomadPop}" \
		--max-af "${gnomadAF}" TMP/jeter1.vcf > TMP/jeter2.vcf

mv TMP/jeter2.vcf TMP/jeter1.vcf


bcftools view TMP/jeter1.vcf -O b -o annot.bcf

bcftools index annot.bcf

###############################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">Annotation in ${bed}</entry>
        <entry key="vcf">${vcf}</entry>
        <entry key="bed">${bed}</entry>
        <entry key="so.acn">${soacn}</entry>
        <entry key="gnomad.pop">${gnomadPop}</entry>
        <entry key="gnomad.af">${gnomadAF}</entry>
        <entry key="gnomad.vcf">${gnomadVcf}</entry>	
        <entry key="versions">${getVersionCmd("bcftools jvarkit/vcffilterso jvarkit/vcfgnomad snpeff")}</entry>
</properties>
EOF
"""
}
