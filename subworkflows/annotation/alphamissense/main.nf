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
include {k1_signature} from '../../modules/utils/k1.nf'


workflow ANNOTATE_ALPHAMISSENSE {
	take:
		fasta
		fai
		dict
		vcfs /* meta, vcf,vcf_index */
	main:
		source_ch = DOWNLOAD(fasta,fai,dict)
		doc_ch = MAKE_DOC(source_ch.bed)
		annotate_ch = ANNOTATE(source_ch.bed, source_ch.tbi,source_ch.header,vcfs)
	emit:
		output = annotate_ch.output
		doc = doc_ch.output
	}


process DOWNLOAD {
tag "${fasta.name}"
afterScript "rm -rf TMP"
label "process_quick"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
        path(fasta)
        path(fai)
        path(dict)
output:
	path("*.tsv.gz"),emit:bed
	path("*.tsv.gz.tbi"),emit:tbi
	path("*.header"),emit:header
script:
    def k1 = k1_signature()
   	def TAG = "ALPHAMISSENSE"
    def whatis="Data from AlphaMissense https://www.science.org/doi/10.1126/science.adg7492"
    def base="https://storage.googleapis.com/dm_alphamissense/AlphaMissense_"
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP

cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg38}    ${base}hg38.tsv.gz
1:${k1.hg19}    ${base}hg19.tsv.gz
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sed 's/^chr//' | sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv
join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv | sort | uniq > TMP/jeter.url

test -s TMP/jeter.url

wget -O TMP/jeter.tsv.gz `cat TMP/jeter.url`

gunzip -c TMP/jeter.tsv.gz  |\\
	grep -v '^#' |\\
	cut -f 1 | uniq | LC_ALL=C sort -T TMP | uniq |\\
	awk '{printf("%s\t%s\\n",\$1,\$1);}' |\\
	jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP bedrenamechr -f "${fasta}" --column 2 --convert SKIP |\\
	awk -F '\t' '{printf("s|^%s\t|%s\t|\\n",\$1,\$2);}' > TMP/jeter.sed


gunzip -c TMP/jeter.tsv.gz |\
	grep -v '^#'  |\
	cut -f 1-4,9,10 |\
	sed -f TMP/jeter.sed |\
	LC_ALL=C sort --buffer-size=1000M -T TMP -t '\t' -k1,1 -k2,2n |\
	uniq |\
	bgzip >  TMP/${TAG}.tsv.gz && \

tabix -s 1 -b 2 -e 2  TMP/${TAG}.tsv.gz

mv TMP/${TAG}.tsv.gz ./
mv TMP/${TAG}.tsv.gz.tbi ./

echo '##INFO=<ID=${TAG}_PATHOGENOCITY,Number=1,Type=Float,Description="${whatis}.">' >  ${TAG}.header
echo '##INFO=<ID=${TAG}_CLASS,Number=.,Type=String,Description="${whatis}.">' >> ${TAG}.header
"""
}


process MAKE_DOC {
executor "local"
input;
    path(bed)
output:
	path("alphamissense.md"),emit:output
script:
    def TAG=bed.name;
"""
cat << EOF > alphamissense.md
Data from AlphaMissense https://www.science.org/doi/10.1126/science.adg7492
EOF
"""
}

process ANNOTATE {
tag "${vcf.name}"
afterScript "rm -rf TMP"
label "process_quick"
conda "${moduleDir}/../../conda/bioinfo.01.yml"

input:
	path(tabix)
	path(tbi)
	path(header)
	tuple val(meta),path(vcf),path(vcf_idx)
output:
    tuple val(meta),path("*.bcf"),path("*.csi"),emit:output
script:
   	def TAG = "ALPHAMISSENSE";
"""
hostname 1>&2
mkdir -p TMP OUTPUT

bcftools annotate \\
    --threads ${task.cpus} \\
    -a "${tabix}" \\
    -h "${header}" \\
    -c "CHROM,POS,REF,ALT,${TAG}_PATHOGENOCITY,${TAG}_CLASS" \\
    --merge-logic "${TAG}_PATHOGENOCITY:max,${TAG}_CLASS:unique" \\
    -O b -o TMP/${TAG}.${vcf.getSimpleName()}.bcf '${vcf}'
 
bcftools \\
    --threads ${task.cpus} \\
    -f index TMP/${TAG}.${vcf.getSimpleName()}.bcf

mv TMP/*.bcf ./
mv TMP/*.bcf.csi ./
"""
}
