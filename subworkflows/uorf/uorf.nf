/*

Copyright (c) 2022 Pierre Lindenbaum

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
include {getKeyValue;getModules} from '../../modules/utils/functions.nf'
include {DOWNLOAD_GTF_01} from '../../modules/gtf/download.gtf.01.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf'
include {GTF_TO_UORF_01} from '../../modules/uorf/gtf.uorf.01.nf'
include {SNPEFF_BUILD_GTF} from '../../modules/snpeff/snpeff.build.gtf.nf'
include {SQRT_FILE} from '../../modules/utils/sqrt.nf'

workflow UORF {
	take:
		meta
		reference
		vcf
		bed
	main:
		version_ch  = Channel.empty()

		gtf_ch = DOWNLOAD_GTF_01(meta, reference )
		version_ch = version_ch.mix(gtf_ch.version)

		vcf2bed_ch = VCF_TO_BED(meta,vcf)
		version_ch = version_ch.mix(vcf2bed_ch.version)


		uorf_ch = GTF_TO_UORF_01(meta, reference, gtf_ch.gtf)
		version_ch = version_ch.mix(uorf_ch.version)

		vep_ch = ANNOTATE_VEP(meta,reference,bed, uorf_ch.gtf,vcf2bed_ch.bed.splitCsv(header: false,sep:'\t',strip:true))

		//SQRT_FILE(meta,COLLECT_TO_FILE(meta,vep_ch.vcf.collect()).output)

		/*
		snpeffdb_ch = SNPEFF_BUILD_GTF(meta, reference, uorf_ch.gtf)
		version_ch = version_ch.mix(snpeffdb_ch.version)

		bed_vcf_snpeff_ch = snpeffdb_ch.snpeffdb.combine(
			vcf2bed_ch.bed.splitCsv(header: false,sep:'\t',strip:true))
		ANNOTATE(meta,uorf_ch.gtf,bed_vcf_snpeff_ch)
		*/
	}

process ANNOTATE_VEP {
tag "${contig}:${start}-${end} ${file(vcf).name}"
memory "5g"
maxForks 1
input:
	val(meta)
	val(reference)
	val(bed)
	val(gtf)
	tuple val(contig),val(start),val(end),val(vcf)
output:
	tuple val(contig),val(start),val(end),path("vep.bcf"),emit:vcf
	path("vep.bcf.csi"),emit:index
script:
"""
hostname 1>&2
set -o pipefail
module load ${getModules("ensembl-vep/104.3 bcftools bedtools")}
mkdir TMP


if [ ! -z "${bed}" ] ; then

	awk -F '\t' '(\$1 == "${contig}"  && !(int(\$2)-10 > ${end} || int(\$3)+10 < ${start} ))' "${bed}" |\
		LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\
		bedtools merge > TMP/user.bed
fi

tabix "${gtf}" "${contig}:${(start as int)+1}-${end}" |\
	awk -F '\t' '(\$3=="gene") {printf("%s\t%d\t%s\\n",\$1,int(\$4)-1,\$5);}' |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n \
	${bed.isEmpty()?"":"bedtools intersect -a - -u -b TMP/user.bed |"} \
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\
	bedtools merge > TMP/jeter.bed

if [ ! -s TMP/jeter.bed ] ; then
	echo "${contig}\t0\t1" > TMP/jeter.bed
fi

bcftools view --regions-file TMP/jeter.bed "${vcf}" |\
        	vep --verbose --format vcf --force_overwrite --gtf "${gtf}"  --fasta "${reference}" --use_given_ref --vcf |\
		bcftools view  -i 'INFO/CSQ ~ "missense"' -O b -o TMP/jeter.bcf

bcftools index TMP/jeter.bcf

mv TMP/jeter.bcf vep.bcf
mv TMP/jeter.bcf.csi vep.bcf.csi
"""
}

process ANNOTATE_SNPEFF {
tag "${contig}:${start}-${end} ${file(vcf).name} / ${dbName}"
memory "5g"
maxForks 1
input:
	val(meta)
	val(gtf)
	tuple val(dbName),val(dbCfg),val(contig),val(start),val(end),val(vcf)
output:
	path("snpeff.bcf"),emit:vcf
	path("snpeff.bcf.csi"),emit:index
script:
"""
hostname 1>&2
set -o pipefail
module load ${getModules("snpEff bcftools bedtools")}
mkdir TMP


tabix "${gtf}" "${contig}:${(start as int)+1}-${end}" |\
	awk -F '\t' '(\$3=="gene") {printf("%s\t%d\t%s\\n",\$1,int(\$4)-1,\$5);}' |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\
	bedtools merge > TMP/jeter.bed

if [ ! -s TMP/jeter.bed ] ; then
	echo "${contig}\t0\t1" > TMP/jeter.bed
fi

bcftools view --regions-file TMP/jeter.bed "${vcf}" |\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${SNPEFF_JAR} eff \
		-config "${dbCfg}" -interval TMP/jeter.bed -nodownload -noNextProt -noMotif -noInteraction -noLog -noStats -chr chr -i vcf -o vcf ${dbName} |\
	bcftools view  -i 'INFO/ANN ~ "missense"' -O b -o TMP/jeter.bcf

bcftools index TMP/jeter.bcf

mv TMP/jeter.bcf snpeff.bcf
mv TMP/jeter.bcf.csi snpeff.bcf.csi
"""
}
