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

def k1_hg19=249250621
def k1_hg38=248956422


workflow ANNOTATE_SV_VCF_01 {
	take:
		reference //bag fasta,fai/dict
		vcf_ch // channel containing [vcf,vcfidx]
	main:
		gtf1_ch = PROCESS_GTF1()
		gtf2_ch = PROCESS_GTF2()
		gnomad_ch = DOWNLOAD_GNOMAD_SV(reference)
		dgv_ch = DOWNLOAD_DGV(reference) 
		ensemblreg_ch = DOWNLOAD_ENSEMBL_REG(reference)

		
		ann_ch= ANNOTATE(
			vcf_ch,
			gtf1_ch.output,
			gtf2_ch.output,
			gnomad_ch.output,
			dgv_ch.output,
			ensemblreg_ch.output
			)

	emit:
		output = ann_ch.output
	}

process PROCESS_GTF1 {
label "process_quick"
output:
	path("gtf.features.*"),emit:output
script:
	def xxxstream=1000
"""
module load htslib
set -xe
set -o pipefail

gunzip -c ${params.gtf} |\\
	awk '/^#/ {next;} {printf("%s\t%d\t%s\t%s\\n",\$1,int(\$4)-1,\$5,\$3); if(\$3=="gene") {  {X=${xxxstream}; printf("%s\t%s\t%s\t%s_%s\\n",\$1,(int(\$4)-X<0?0:int(\$4)-X),\$4,(\$7=="+"?"upstream":"downstream"),X); printf("%s\t%s\t%s\t%s_%s\\n",\$1,\$5,\$5+X,(\$7=="+"?"downstream":"upstream"),X);}  } }' |\\
	LC_ALL=C sort -T . -t '\t' -k1,1 -k2,2n |\\
	uniq |\\
	bgzip > gtf.features.bed.gz

tabix -p bed -f gtf.features.bed.gz

echo '##INFO=<ID=GTF_FEATURE,Number=.,Type=String,Description="features from ${params.gtf}">' > gtf.features.hdr
"""
}



process PROCESS_GTF2 {
label "process_quick"
output:
	path("gtf.genes.*"),emit:output
script:
"""
set -xe
module load htslib jvarkit
set -o pipefail


java -jar \${JVARKIT_DIST}/jvarkit.jar gtf2bed --columns "gtf.feature,gene_name" "${params.gtf}" |\\
	awk -F '\t' '(\$4=="gene" && \$5!=".")' |\\
	cut -f1,2,3,5 |\\
	LC_ALL=C sort -T . -t '\t' -k1,1 -k2,2n |\\
	uniq | bgzip > gtf.genes.bed.gz

tabix -p bed -f gtf.genes.bed.gz

echo '##INFO=<ID=GENE,Number=.,Type=String,Description="gene from ${params.gtf}">' > gtf.genes.hdr

"""
}

process DOWNLOAD_GNOMAD_SV {
label "process_quick"
input:
	path(reference)
output:
	path("gnomad.sv.*"),emit:output
script:
	def fai = reference.find{it.name.endsWith(".fai")}.first()
	def dict = reference.find{it.name.endsWith(".dict")}.first()
"""
module load htslib jvarkit
set -xe

mkdir -p TMP
cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1_hg38}	https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.bed.gz
1:${k1_hg19}	https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.bed.gz
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sed 's/^chr//' | sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv

join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv | sort | uniq > TMP/jeter.url

test -s TMP/jeter.url

wget -O - `cat TMP/jeter.url` |\\
	gunzip -c |\\
	awk -F '\t' '(NR==1){col=-1;for(i=1;i<=NF;i++) {if(col<0 && (\$i=="GRPMAX_AF" ||\$i=="POPMAX_AF")) col=i;} next;} {if(col>1) printf("%s\\t%s\\t%s\\t%s\\t%s\\n",\$1,\$2,\$3,\$4,(\$col ~ /^[0-9Ee\\-\\.]+\$/? \$col:"."));}'  |\
	java -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${dict}" --column 1 --convert SKIP  |\\
	LC_ALL=C sort -T . -t '\t' -k1,1 -k2,2n |\\
        uniq | bgzip > gnomad.sv.bed.gz
	
tabix -f -p bed gnomad.sv.bed.gz


echo '##INFO=<ID=GNOMAD_ID,Number=1,Type=String,Description="First GNOMAD SV ID found for this variant">' >  gnomad.sv.hdr
echo '##INFO=<ID=GNOMAD_AF,Number=1,Type=Float,Description="GNOMAD SV MAX AF">' >> gnomad.sv.hdr

"""	
}

process DOWNLOAD_ENSEMBL_REG {
label "process_quick"
input:
	path(reference)
output:
	path("ensembl.reg.*"),emit:output
script:
	def fai = reference.find{it.name.endsWith(".fai")}.first()
	def dict = reference.find{it.name.endsWith(".dict")}.first()
"""
hostname 1>&2
set -xe
module load htslib jvarkit

mkdir -p TMP
cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1_hg38}	https://ftp.ensembl.org/pub/release-111/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz
1:${k1_hg19}	https://ftp.ensembl.org/pub/grch37/current/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.gff.gz
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sed 's/^chr//' | sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv

join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv | sort | uniq > TMP/jeter.url

test -s TMP/jeter.url

wget -O - `cat TMP/jeter.url` |\\
	gunzip -c |\
	awk '/^#/ {next;} {printf("%s\t%d\t%s\t%s\\n",\$1,int(\$4)-1,\$5,\$3);}' |\\
	java -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${dict}" --column 1 --convert SKIP  |\\
	LC_ALL=C sort -T . -t '\t' -k1,1 -k2,2n | uniq |\\
	bgzip > ensembl.reg.bed.gz

tabix -p bed -f ensembl.reg.bed.gz

echo '##INFO=<ID=ENS_REG,Number=.,Type=String,Description="features from Ensembl Regulation">' > ensembl.reg.hdr

"""
}

process DOWNLOAD_DGV {
label "process_quick"
input:
	path(reference)
output:
       	path("dgv.*"),emit:output
script:
	def fai = reference.find{it.name.endsWith(".fai")}.first()
	def dict = reference.find{it.name.endsWith(".dict")}.first()
"""
hostname 1>&2
module load htslib jvarkit
set -xe

mkdir -p TMP
cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1_hg38}	http://dgv.tcag.ca/dgv/docs/GRCh38_hg38_variants_2020-02-25.txt
1:${k1_hg19}	http://dgv.tcag.ca/dgv/docs/GRCh37_hg19_variants_2020-02-25.txt
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sed 's/^chr//' | sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv

join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv | sort | uniq > TMP/jeter.url

test -s TMP/jeter.url

wget  -O - `cat TMP/jeter.url` |\
	awk -F '\t' '(NR==1){col=-1;for(i=1;i<=NF;i++) {if(col<0 && \$i=="frequency") col=i;} next;} {if(col>1) printf("%s\\t%s\\t%s\\t%s\\t%s\\n",\$2,\$3,\$4,\$1,(\$col ~ /^[0-9Ee\\-\\.]+\$/? \$col:"."));}'  |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${dict}" --column 1 --convert SKIP |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\\
	bgzip > dgv.bed.gz

tabix  -p bed -f dgv.bed.gz


echo '##INFO=<ID=DGV_ID,Number=1,Type=String,Description="First DGV ID found for the variant">' >  dgv.hdr
echo '##INFO=<ID=DGV_AF,Number=1,Type=Float,Description="DGV_FREQUENCY">' >> dgv.hdr
"""
}



process ANNOTATE {
label "process_quick"
tag "${vcf.name}"
afterScript "rm -rf TMP"
input:
	tuple path(vcf),path(vcfidx)
	path(gtf1_files)
	path(gtf2_files)
	path(gnomad_files)
	path(dgv_files)
	path(ensemblreg_files)
output:
	tuple path("${vcf.getBaseName()}.ann.bcf"),path("${vcf.getBaseName()}.ann.bcf.csi"),emit:output
script:
	def overlap=0.7
"""
hostname 1>&2
module load bcftools
mkdir -p TMP
set -x
                
bcftools annotate --force --annotations "${gtf1_files.find{it.name.endsWith(".gz")}.first()}" \
	-h "${gtf1_files.find{it.name.endsWith(".hdr")}.first()}" \
	-c "CHROM,FROM,TO,GTF_FEATURE" \
	--merge-logic GTF_FEATURE:unique \
	-O u -o TMP/jeter2.bcf "${vcf}"
                
mv TMP/jeter2.bcf TMP/jeter1.bcf

bcftools annotate --force --annotations "${gtf2_files.find{it.name.endsWith(".gz")}.first()}" \
	-h "${gtf2_files.find{it.name.endsWith(".hdr")}.first()}" \
	-c "CHROM,FROM,TO,GENE" \
	--merge-logic GENE:unique \
	-O u -o TMP/jeter2.bcf TMP/jeter1.bcf
                
mv TMP/jeter2.bcf TMP/jeter1.bcf

bcftools annotate --force --annotations "${ensemblreg_files.find{it.name.endsWith(".gz")}.first()}" \
	-h "${ensemblreg_files.find{it.name.endsWith(".hdr")}.first()}" \
	-c "CHROM,FROM,TO,ENS_REG" \
	--merge-logic ENS_REG:unique \
	-O u -o TMP/jeter2.bcf TMP/jeter1.bcf
                
mv TMP/jeter2.bcf TMP/jeter1.bcf


bcftools annotate --force --annotations "${gnomad_files.find{it.name.endsWith(".gz")}.first()}" \
	-h "${gnomad_files.find{it.name.endsWith(".hdr")}.first()}" \\
	-c "CHROM,FROM,TO,GNOMAD_ID,GNOMAD_AF" \\
	--merge-logic 'GNOMAD_ID:first,GNOMAD_AF:max' \\
	--min-overlap "${overlap}:${overlap}" \\
	-O u -o TMP/jeter2.bcf TMP/jeter1.bcf
                
mv TMP/jeter2.bcf TMP/jeter1.bcf

bcftools annotate --force --annotations "${dgv_files.find{it.name.endsWith(".gz")}.first()}" \
	-h "${dgv_files.find{it.name.endsWith(".hdr")}.first()}" \\
	-c "CHROM,FROM,TO,DGV_ID,DGV_AF" \\
	--merge-logic 'DGV_ID:first,DGV_AF:max' \\
	--min-overlap "${overlap}:${overlap}" \\
	-O u -o TMP/jeter2.bcf TMP/jeter1.bcf
                
mv TMP/jeter2.bcf TMP/jeter1.bcf




bcftools index -f TMP/jeter1.bcf

mv TMP/jeter1.bcf ${vcf.getBaseName()}.ann.bcf
mv TMP/jeter1.bcf.csi ${vcf.getBaseName()}.ann.bcf.csi

"""
}
