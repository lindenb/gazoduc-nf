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

workflow ANNOTATE_SV {
	take:
		meta 
		fasta
		fai
		dict
		gtf
		vcf_ch // channel containing [vcf,vcfidx]
	main:
        versions = Channel.empty()
		PROCESS_GTF1(gtf)
		PROCESS_GTF2(gtf)
		DOWNLOAD_GNOMAD_SV(dict)
		versions = versions.mix(DOWNLOAD_GNOMAD_SV.out.versions)

		DGV_DOWNLOAD(dict)
		versions = versions.mix(DGV_DOWNLOAD.out.versions)
		
		ENSEMBL_REG_DOWNLOAD(dict)
		versions = versions.mix(ENSEMBL_REG_DOWNLOAD.out.versions)

		DOWNLOAD_DECODE_LRS(dict)

		
		ann_ch= ANNOTATE(
			vcf_ch,
			PROCESS_GTF1.out.bed,
			PROCESS_GTF2.out.bed,
			DOWNLOAD_GNOMAD_SV.out.vcf,
			DGV_DOWNLOAD.out.bed,
			ENSEMBL_REG_DOWNLOAD.out.bed,
			DOWNLOAD_DECODE_LRS.out.tabix
			)

		ch1 = fasta.branch{v->
			fasta_hg38: v[0].ucsc_name!=null && v[0].ucsc_name=="hg38"
			others: true
			}

		vcf = Channel.empty()
			
		VEP_GRCH38(
			meta,
			ch1.fasta_hg38,
			fai,
			dict,
			vcf_ch
			)
		vcf = vcf.mix(VEP_GRCH38.out.vcf)

	emit:
		vcf
        versions
	}

process PROCESS_GTF1 {
tag "${gtf.name}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(gtf),path(gtf_tbi)
output:
	tuple val(meta1),path("*.bed.gz"),path("*.bed.gz.tbi"),path("*.hdr"),emit:bed
script:
	def xxxstream=task.ext.extend?:1000
"""
set -xe
set -o pipefail

gunzip -c ${params.gtf} |\\
	awk '/^#/ {next;} {printf("%s\t%d\t%s\t%s\\n",\$1,int(\$4)-1,\$5,\$3); if(\$3=="gene") {  {X=${xxxstream}; printf("%s\t%s\t%s\t%s_%s\\n",\$1,(int(\$4)-X<0?0:int(\$4)-X),\$4,(\$7=="+"?"upstream":"downstream"),X); printf("%s\t%s\t%s\t%s_%s\\n",\$1,\$5,\$5+X,(\$7=="+"?"downstream":"upstream"),X);}  } }' |\\
	LC_ALL=C sort -T . -t '\t' -k1,1 -k2,2n |\\
	uniq |\\
	bgzip > gtf.features.bed.gz

tabix -p bed -f gtf.features.bed.gz

echo '##INFO=<ID=GTF_FEATURE,Number=.,Type=String,Description="features from ${gtf}">' > gtf.features.hdr
"""
}



process PROCESS_GTF2 {
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(gtf),path(gtf_tbi)
output:
	tuple val(meta1),path("*.bed.gz"),path("*.bed.gz.tbi"),path("*.hdr"),emit:bed
script:
"""
set -xe
set -o pipefail


jvarkit gtf2bed --columns "gtf.feature,gene_name" "${gtf}" |\\
	awk -F '\t' '(\$4=="gene" && \$5!=".")' |\\
	cut -f1,2,3,5 |\\
	LC_ALL=C sort -T . -t '\t' -k1,1 -k2,2n |\\
	uniq | bgzip > gtf.genes.bed.gz

tabix -p bed -f gtf.genes.bed.gz

echo '##INFO=<ID=GENE,Number=.,Type=String,Description="gene from ${gtf}">' > gtf.genes.hdr

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
script:
	def url = "";
	if(meta.ucsc_name=="hg38") {
		url = "https://github.com/DecodeGenetics/LRS_SV_sets/raw/refs/heads/master/ont_sv_high_confidence_SVs.sorted.vcf.gz"
		}
	
"""
mkdir -p TMP


if ${!url.isEmpty()}
then
	curl -L "${url}" |\\
		bcftools query -f '%CHROM\t%POS0\t%END\t%ID\\_%SVTYPE\\n'  |\
		jvarkit bedrenamechr -f "${dict}" --column 1 --convert SKIP  |\\
		LC_ALL=C sort -T . -t '\t' -k1,1 -k2,2n | uniq |\\
		bgzip >  decode.bed.gz

else
	touch decode.bed
	bgzip decode.bed
fi


tabix -p bed -f decode.bed.gz

echo '##INFO=<ID=DECODE_LRS,Number=.,Type=String,Description="DECODE SV. ${url} Accompanies Long-read sequencing of 3,622 Icelanders provides insight into the role of structural variants in human diseases and other traits Beyter, D. et al., Nature Genetics, 2021">' > decode.hdr
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
	tuple val(meta3),path(gnomad_vcf),path(gnomad_tbi)
	tuple val(meta4),path(dgv_bed),path(dgv_tbi),path(dgv_hdr)
	tuple val(meta5),path(reg_bed),path(reg_tbi),path(reg_hdr)
	tuple val(meta6),path(decode_bed),path(decode_tbi),path(decode_hdr)
output:
	tuple val(meta),path("*.bcf"),path("*.bcf.csi"),emit:vcf
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
	--force --annotations "${gnomad_bed}" \
	-h "${gnomad_hdr}" \\
	-c "CHROM,FROM,TO,GNOMAD_ID,GNOMAD_AF" \\
	--merge-logic 'GNOMAD_ID:first,GNOMAD_AF:max' \\
	--min-overlap "${overlap}:${overlap}" \\
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


bcftools annotate  --threads ${task.cpus} \\
	--force --annotations "${dgv_bed}" \
	-h "${dgv_hdr}" \\
	-c "CHROM,FROM,TO,DGV_ID,DGV_AF" \\
	--merge-logic 'DGV_ID:first,DGV_AF:max' \\
	--min-overlap "${overlap}:${overlap}" \\
	-O b9 -o TMP/jeter2.bcf TMP/jeter1.bcf
                
mv TMP/jeter2.bcf TMP/jeter1.bcf


bcftools index -f TMP/jeter1.bcf

mv TMP/jeter1.bcf ${prefix}.bcf
mv TMP/jeter1.bcf.csi ${prefix}.bcf.csi

"""
}
