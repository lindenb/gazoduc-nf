/*

Copyright (c) 2026 Pierre Lindenbaum

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
include { VEP as VEP_GRCH38       } from './grch38.nf'
include { DOWNLOAD_GNOMAD_SV      } from '../../../modules/gnomad_sv/download.vcf'
include { DGV_DOWNLOAD            } from '../../../modules/dgv/download'
include { ENSEMBL_REG_DOWNLOAD    } from '../../../modules/ensemblreg/download'
include { JVARKIT_VCFGNOMADSV     } from '../../../modules/jvarkit/vcfgnomadsv'

workflow ANNOTATE_SV {
	take:
		metadata
		fasta
		fai
		dict
		gtf
		vcf_ch // channel containing [vcf,vcfidx]
	main:
        versions = Channel.empty()
		multiqc  = Channel.empty()
		
		PROCESS_GTF1(gtf)
		versions = versions.mix(PROCESS_GTF1.out.versions)
		PROCESS_GTF2(gtf)
		versions = versions.mix(PROCESS_GTF2.out.versions)
		
		DOWNLOAD_GNOMAD_SV(dict)
		versions = versions.mix(DOWNLOAD_GNOMAD_SV.out.versions)
		
		
		JVARKIT_VCFGNOMADSV(
			DOWNLOAD_GNOMAD_SV.out.vcf,
			vcf_ch.map{meta,vcf,tbi->[meta,vcf]}
			)
		versions = versions.mix(JVARKIT_VCFGNOMADSV.out.versions)
		vcf_ch = JVARKIT_VCFGNOMADSV.out.vcf

		DGV_DOWNLOAD(dict)
		versions = versions.mix(DGV_DOWNLOAD.out.versions)
		
		ENSEMBL_REG_DOWNLOAD(dict)
		versions = versions.mix(ENSEMBL_REG_DOWNLOAD.out.versions)

		DOWNLOAD_DECODE_LRS(dict)
		versions = versions.mix(DOWNLOAD_DECODE_LRS.out.versions)
		
		ANNOTATE(
			vcf_ch,
			PROCESS_GTF1.out.bed,
			PROCESS_GTF2.out.bed,
			DGV_DOWNLOAD.out.bed,
			ENSEMBL_REG_DOWNLOAD.out.bed,
			DOWNLOAD_DECODE_LRS.out.tabix
			)
		versions = versions.mix(ANNOTATE.out.versions)
		vcf_ch = ANNOTATE.out.vcf

		ch1 = fasta.branch{v->
			fasta_hg38: v[0].ucsc_name!=null && v[0].ucsc_name=="hg38"
			others: true
			}
		
		vcf = Channel.empty()
			
		VEP_GRCH38(
			metadata,
			ch1.fasta_hg38,
			fai,
			dict,
			vcf_ch
			)
		vcf = vcf.mix(VEP_GRCH38.out.vcf)

	emit:
		vcf
        versions
		multiqc
	}

process PROCESS_GTF1 {
tag "${meta.id}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta),path(gtf)
output:
	tuple val(meta),path("*.bed.gz"),path("*.bed.gz.tbi"),path("*.hdr"),emit:bed
	path("versions.yml"),emit: versions
script:
	def xxxstream= task.ext.extend?:1000
	def prefix= task.ext.prefix?:"${meta.id}.features"
"""
gunzip -c "${gtf}"|\\
	awk '/^#/ {next;} {printf("%s\t%d\t%s\t%s\\n",\$1,int(\$4)-1,\$5,\$3); if(\$3=="gene") {  {X=${xxxstream}; printf("%s\t%s\t%s\t%s_%s\\n",\$1,(int(\$4)-X<0?0:int(\$4)-X),\$4,(\$7=="+"?"upstream":"downstream"),X); printf("%s\t%s\t%s\t%s_%s\\n",\$1,\$5,\$5+X,(\$7=="+"?"downstream":"upstream"),X);}  } }' |\\
	LC_ALL=C sort -T . -t '\t' -k1,1 -k2,2n |\\
	uniq |\\
	bgzip > ${prefix}.bed.gz

tabix -p bed -f ${prefix}.bed.gz

echo '##INFO=<ID=GTF_FEATURE,Number=.,Type=String,Description="features from ${gtf}">' > ${prefix}.hdr

touch versions.yml
"""
stub:
	def prefix= task.ext.prefix?:"${meta.id}.features"
"""
touch  versions.yml ${prefix}.bed.gz ${prefix}.bed.gz.tbi ${prefix}.hdr
"""
}



process PROCESS_GTF2 {
tag "${meta.id}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta),path(gtf)
output:
	tuple val(meta),path("*.bed.gz"),path("*.bed.gz.tbi"),path("*.hdr"),emit:bed
	path("versions.yml"),emit: versions
script:
	def prefix= task.ext.prefix?:"${meta.id}.genes"
"""
set -xe
set -o pipefail


jvarkit gtf2bed --columns "gtf.feature,gene_name,gene_id" "${gtf}" |\\
	awk -F '\t' '(\$4=="gene" && \$5!="." ) {OFS="\t";if(\$5=="." || $\5=="") \$5=\$6; if(\$5=="." || $\5=="") next; print;}' |\\
	cut -f1,2,3,5 |\\
	LC_ALL=C sort -T . -t '\t' -k1,1 -k2,2n |\\
	uniq | bgzip > ${prefix}.bed.gz

tabix -p bed -f ${prefix}.bed.gz

echo '##INFO=<ID=GENE,Number=.,Type=String,Description="gene from ${gtf}">' > ${prefix}.hdr

touch versions.yml
"""
stub:
	def prefix= task.ext.prefix?:"${meta.id}.genes"
"""
touch  versions.yml ${prefix}.bed.gz ${prefix}.bed.gz.tbi ${prefix}.hdr
"""
}

process DOWNLOAD_DECODE_LRS {
tag "${meta.id?:""}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta),path(dict)
output:
	tuple val(meta),path("*.bed.gz"),path("*.bed.gz.tbi"),path("*.hdr"),emit:tabix
	path("versions.yml"),emit: versions
script:
	def url = "";
	if(meta.ucsc_name=="hg38") {
		url = "https://github.com/DecodeGenetics/LRS_SV_sets/raw/refs/heads/master/ont_sv_high_confidence_SVs.sorted.vcf.gz"
		}
	def prefix = task.ext.prefix?:"decode.lrs"
"""
mkdir -p TMP


if ${!url.isEmpty()}
then
	curl -L "${url}" |\\
		bcftools query -f '%CHROM\t%POS0\t%END\t%ID\\_%SVTYPE\\n'  |\
		jvarkit bedrenamechr -f "${dict}" --column 1 --convert SKIP  |\\
		LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n | uniq |\\
		bgzip >  ${prefix}.bed.gz

else
	touch ${prefix}.bed
	bgzip ${prefix}.bed
fi


tabix -p bed -f ${prefix}.bed.gz

echo '##INFO=<ID=DECODE_LRS,Number=.,Type=String,Description="DECODE SV. ${url} Accompanies Long-read sequencing of 3,622 Icelanders provides insight into the role of structural variants in human diseases and other traits Beyter, D. et al., Nature Genetics, 2021">' > ${prefix}.hdr
touch versions.yml
"""
stub:
	def prefix = task.ext.prefix?:"decode.lrs"
"""
touch  versions.yml ${prefix}.bed.gz ${prefix}.bed.gz.tbi ${prefix}.hdr
"""
}




process ANNOTATE {
tag "${meta.id}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta ),path(vcf),path(vcfidx)
	tuple val(meta1),path(gtf1_bed),path(gtf1_tbi),path(gtf1_hdr)
	tuple val(meta2),path(gtf2_bed),path(gtf2_tbi),path(gtf2_hdr)
	tuple val(meta4),path(dgv_bed),path(dgv_tbi)
	tuple val(meta5),path(reg_bed),path(reg_tbi),path(reg_hdr)
	tuple val(meta6),path(decode_bed),path(decode_tbi),path(decode_hdr)
output:
	tuple val(meta),path("*.vcf"),path("*.vcf.gz"),emit:vcf
	path("versions.yml"),emit:versions
script:
	def overlap= task.ext.overlap?:0.7
	def prefix = task.ext.prefix?:vcf.baseName+".ann"
"""
hostname 1>&2
mkdir -p TMP
set -x
                
bcftools annotate --threads ${task.cpus} \\
	--force --annotations "${gtf1_bed}" \\
	-h "${gtf1_hdr}" \\
	-c "CHROM,FROM,TO,GTF_FEATURE" \\
	--merge-logic GTF_FEATURE:unique \\
	-O u -o TMP/jeter2.bcf "${vcf}"
                
mv TMP/jeter2.bcf TMP/jeter1.bcf

bcftools annotate  --threads ${task.cpus} \\
	--force --annotations "${gtf2_bed}" \\
	-h "${gtf2_hdr}" \\
	-c "CHROM,FROM,TO,GENE" \\
	--merge-logic GENE:unique \\
	-O u -o TMP/jeter2.bcf TMP/jeter1.bcf
                
mv TMP/jeter2.bcf TMP/jeter1.bcf

bcftools annotate  --threads ${task.cpus} \\
	--force --annotations "${reg_bed}" \
	-h "${reg_hdr}" \\
	-c "CHROM,FROM,TO,ENS_REG" \\
	--merge-logic ENS_REG:unique \\
	-O u -o TMP/jeter2.bcf TMP/jeter1.bcf
                
mv TMP/jeter2.bcf TMP/jeter1.bcf


bcftools annotate  --threads ${task.cpus} \\
	--force --annotations "${decode_bed}" \
	-h "${decode_hdr}" \\
	-c "CHROM,FROM,TO,DECODE_LRS" \\
	--merge-logic 'DECODE_LRS:unique' \\
	--min-overlap "${overlap}:${overlap}" \\
	-O u -o TMP/jeter2.bcf TMP/jeter1.bcf
                
mv TMP/jeter2.bcf TMP/jeter1.bcf


cat << EOF  > TMP/dgv.hdr
##INFO=<ID=DGV_ID,Number=1,Type=String,Description="ID in DGV">
##INFO=<ID=DGV_AF,Number=1,Type=Float,Description="Freq in DGV">
EOF

bcftools annotate  --threads ${task.cpus} \\
	--force --annotations "${dgv_bed}" \
	-h TMP/dgv.hdr \\
	-c "CHROM,FROM,TO,DGV_ID,DGV_AF" \\
	--merge-logic 'DGV_ID:first,DGV_AF:max' \\
	--min-overlap "${overlap}:${overlap}" \\
	-O b9 -o TMP/jeter2.bcf TMP/jeter1.bcf
                
mv TMP/jeter2.bcf TMP/jeter1.bcf

bcftools view --threads ${task.cpus}  -O z -o TMP/jeter2.vcf.gz TMP/jeter1.bcf
mv  TMP/jeter2.vcf.gz TMP/jeter1.vcf.gz


bcftools index -t -f TMP/jeter1.vcf.gz

mv TMP/jeter1.vcf.gz ${prefix}.vcf.gz
mv TMP/jeter1.vcf.gz.tbi ${prefix}.vcf.gz.tbi

touch versions.yml
"""
}
