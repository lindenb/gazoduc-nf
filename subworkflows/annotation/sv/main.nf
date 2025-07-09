/*

Copyright (c) 2025 Pierre Lindenbaum

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
include {k1_signature         } from '../../../modules/utils/k1.nf'
include {isGRCh38             } from '../../../modules/utils/k1.nf'
include {VEP as VEP_GRCH38    } from './grch38.nf'

def k1 = k1_signature()


workflow ANNOTATE_SV {
	take:
		meta 
		fasta
		fai
		dict
		gtf
		vcf_ch // channel containing [vcf,vcfidx]
	main:
		gtf1_ch = PROCESS_GTF1(gtf)
		gtf2_ch = PROCESS_GTF2(gtf)
		gnomad_ch = DOWNLOAD_GNOMAD_SV(fasta,fai,dict)
		dgv_ch = DOWNLOAD_DGV(fasta,fai,dict) 
		ensemblreg_ch = DOWNLOAD_ENSEMBL_REG(fasta,fai,dict)
		decode_ch = DOWNLOAD_DECODE_LRS(fasta,fai,dict)

		
		ann_ch= ANNOTATE(
			vcf_ch,
			gtf1_ch.output,
			gtf2_ch.output,
			gnomad_ch.output,
			dgv_ch.output,
			ensemblreg_ch.output,
			decode_ch.output
			)
		vcf = ANNOTATE.out.vcf
		if(isGRCh38(fai[1])) {
				VEP_GRCH38(meta,fasta,fai,dict,vcf)
				vcf = VEP_GRCH38.out.vcf
			} else {
				throw new IllegalArgumentException("VEP: unknown build");
			}
	emit:
		vcf
	}

process PROCESS_GTF1 {
tag "${gtf.name}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(gtf),path(gtf_tbi)
output:
	tuple val(meta1),path("*.bed.gz"),path("*.bed.gz.tbi"),path("*.hdr"),emit:output
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
	tuple val(meta1),path("*.bed.gz"),path("*.bed.gz.tbi"),path("*.hdr"),emit:output
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

process DOWNLOAD_GNOMAD_SV {
tag "${fasta.name}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
label "process_single"
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
output:
	tuple val(meta1),path("*.bed.gz"),path("*.bed.gz.tbi"),path("*.hdr"),emit:output
script:
	def base="https://storage.googleapis.com/gcp-public-data--gnomad"
"""
set -xe

mkdir -p TMP
cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg38}\t${base}/release/4.1/genome_sv/gnomad.v4.1.sv.sites.bed.gz
1:${k1.hg19}\t${base}/papers/2019-sv/gnomad_v2.1_sv.sites.bed.gz
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' |\\
	sed 's/^chr//' |\\
	sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv

join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv | sort | uniq > TMP/jeter.url

test -s TMP/jeter.url

wget -O - `cat TMP/jeter.url` |\\
	gunzip -c |\\
	awk -F '\t' '(NR==1){col=-1;for(i=1;i<=NF;i++) {if(col<0 && (\$i=="GRPMAX_AF" ||\$i=="POPMAX_AF")) col=i;} next;} {if(col>1) printf("%s\\t%s\\t%s\\t%s\\t%s\\n",\$1,\$2,\$3,\$4,(\$col ~ /^[0-9Ee\\-\\.]+\$/? \$col:"."));}'  |\
	jvarkit bedrenamechr -f "${dict}" --column 1 --convert SKIP  |\\
	LC_ALL=C sort -T . -t '\t' -k1,1 -k2,2n |\\
        uniq | bgzip > gnomad.sv.bed.gz
	
tabix -f -p bed gnomad.sv.bed.gz


echo '##INFO=<ID=GNOMAD_ID,Number=1,Type=String,Description="First GNOMAD SV ID found for this variant">' >  gnomad.sv.hdr
echo '##INFO=<ID=GNOMAD_AF,Number=1,Type=Float,Description="GNOMAD SV MAX AF">' >> gnomad.sv.hdr

"""	
}

process DOWNLOAD_DECODE_LRS {
tag "${meta1.id?:""}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
output:
	tuple val(meta1),path("*.bed.gz"),path("*.bed.gz.tbi"),path("*.hdr"),emit:output
script:

"""

mkdir -p TMP
cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg38}\thttps://github.com/DecodeGenetics/LRS_SV_sets/raw/refs/heads/master/ont_sv_high_confidence_SVs.sorted.vcf.gz
EOF


awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sed 's/^chr//' | sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv

join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv | sort | uniq > TMP/jeter.url

if -s TMP/jeter.url
then
	wget -O - `cat TMP/jeter.url` |\\
		bcftools query -f '%CHROM\t%POS0\t%END\t%ID\\_%SVTYPE\\n'  |\
		jvarkit bedrenamechr -f "${dict}" --column 1 --convert SKIP  |\\
		LC_ALL=C sort -T . -t '\t' -k1,1 -k2,2n | uniq |\\
		bgzip >  decode.bed.gz

else
	touch decode.bed
	bgzip decode.bed
fi


tabix -p bed -f decode.bed.gz

echo '##INFO=<ID=DECODE_LRS,Number=.,Type=String,Description="DECODE SV. https://github.com/DecodeGenetics/LRS_SV_sets Accompanies Long-read sequencing of 3,622 Icelanders provides insight into the role of structural variants in human diseases and other traits Beyter, D. et al., Nature Genetics, 2021">' > decode.hdr
"""
}


process DOWNLOAD_ENSEMBL_REG {
tag "${meta1.id?:""}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
output:
	tuple val(meta1),path("*.bed.gz"),path("*.bed.gz.tbi"),path("*.hdr"),emit:output
script:
	def base="https://ftp.ensembl.org/pub"
"""
hostname 1>&2
set -xe

mkdir -p TMP
cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg38}\t${base}/release-111/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz
1:${k1.hg19}\t${base}/grch37/current/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.gff.gz
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sed 's/^chr//' | sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv

join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv | sort | uniq > TMP/jeter.url

test -s TMP/jeter.url

wget -O - `cat TMP/jeter.url` |\\
	gunzip -c |\
	awk '/^#/ {next;} {printf("%s\t%d\t%s\t%s\\n",\$1,int(\$4)-1,\$5,\$3);}' |\\
	jvarkit bedrenamechr -f "${dict}" --column 1 --convert SKIP  |\\
	LC_ALL=C sort -T . -t '\t' -k1,1 -k2,2n | uniq |\\
	bgzip > ensembl.reg.bed.gz

tabix -p bed -f ensembl.reg.bed.gz

echo '##INFO=<ID=ENS_REG,Number=.,Type=String,Description="features from Ensembl Regulation">' > ensembl.reg.hdr
"""
}

process DOWNLOAD_DGV {
tag "${meta1.id?:""}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
output:
	tuple val(meta1),path("*.bed.gz"),path("*.bed.gz.tbi"),path("*.hdr"),emit:output
script:
	def base="http://dgv.tcag.ca/dgv/docs";
"""
hostname 1>&2
set -xe

mkdir -p TMP
cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg38}\t${base}/GRCh38_hg38_variants_2020-02-25.txt
1:${k1.hg19}\t${base}/GRCh37_hg19_variants_2020-02-25.txt
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sed 's/^chr//' | sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv

join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv | sort | uniq > TMP/jeter.url

test -s TMP/jeter.url

wget  -O - `cat TMP/jeter.url` |\
	awk -F '\t' '(NR==1){col=-1;for(i=1;i<=NF;i++) {if(col<0 && \$i=="frequency") col=i;} next;} {if(col>1) printf("%s\\t%s\\t%s\\t%s\\t%s\\n",\$2,\$3,\$4,\$1,(\$col ~ /^[0-9Ee\\-\\.]+\$/? \$col:"."));}'  |\
	java bedrenamechr -f "${dict}" --column 1 --convert SKIP |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\\
	bgzip > dgv.bed.gz

tabix  -p bed -f dgv.bed.gz


echo '##INFO=<ID=DGV_ID,Number=1,Type=String,Description="First DGV ID found for the variant">' >  dgv.hdr
echo '##INFO=<ID=DGV_AF,Number=1,Type=Float,Description="DGV_FREQUENCY">' >> dgv.hdr
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
	tuple val(meta3),path(gnomad_bed),path(gnomad_tbi),path(gnomad_hdr)
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
