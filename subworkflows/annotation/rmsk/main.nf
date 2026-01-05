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
include {k1_signature} from '../../../modules/utils/k1.nf'

workflow RMSK {
	take:
        meta
		fasta
		fai
		dict
		vcfs /* meta, vcf,vcf_index */
	main:
        versions = Channel.empty()
		DOWNLOAD(fasta,fai,dict)
        versions = versions.mix(DOWNLOAD.out.versions)

		ANNOTATE(DOWNLOAD.out.output,vcfs)
        versions = versions.mix(ANNOTATE.out.versions)
	emit:
		vcf = ANNOTATE.out.vcf
		versions
		doc = DOWNLOAD.out.doc
	}


process DOWNLOAD{
tag "${meta1.id?:fasta.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
output:
	tuple val(meta1),path("*.bed.gz"), path("*.bed.gz.tbi"), path("*.header"),emit:output
	path("versions.yml"),emit:versions
	path("doc.md"),emit:doc
script:
	def k1 = k1_signature()
    def base = "https://hgdownload.cse.ucsc.edu/goldenPath"
	def TAG = task.ext.tag?:"RMSK"

"""
hostname 1>&2
mkdir -p TMP



cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg38}\t${base}/hg38/database/rmsk.txt.gz
1:${k1.hg19}\t${base}/hg19/database/rmsk.txt.gz
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' |\\
    sed 's/^chr//' |\\
    sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv
join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv |\\
sort | uniq > TMP/jeter.url

URL=`cat TMP/jeter.url`

set -o pipefail
wget --no-check-certificate -O - "\${URL}" |\\
	gunzip -c |\\
	cut -f6-8 |\\
	jvarkit bedrenamechr -XX:-UsePerfData  -Djava.io.tmpdir=TMP -f "${fasta}" --column 1 --convert SKIP  |\\
		sort  -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n |\\
		bedtools merge |\\
		sed 's/\$/\t1/' |\\
		bgzip > TMP/${TAG}.bed.gz
	

tabix -p bed -f TMP/${TAG}.bed.gz

mv TMP/${TAG}.bed.gz ./
mv TMP/${TAG}.bed.gz.tbi ./

echo '##INFO=<ID=${TAG},Number=0,Type=Flag,Description="Repeat Masker from UCSC">' > ${TAG}.header

cat << 'EOF' > doc.md
# annotations:repeatmasker

`INFO/${TAG}` : Repeat Masker regions from UCSC

EOF

cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "\${URL}"
END_VERSIONS
"""
}



process ANNOTATE {
label "process_single"
tag "${meta.id?:vcf.name}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(tabix),path(tbi),path(header)
	tuple val(meta),path(vcf),path(vcf_idx)
output:
	tuple val(meta), path("*.bcf"), path("*.bcf.csi"),emit:vcf
	path("versions.yml"),emit:versions
script:
	def TAG = task.ext.tag?:"RMSK"
	def prefix=task.ext.prefix?:vcf.baseName+".rmsk"
"""
mkdir -p TMP

bcftools annotate \\
	--threads ${task.cpus} \\
	--write-index \\
	-a "${tabix}" \\
	-h "${header}" \\
	-c "CHROM,FROM,TO,${TAG}"  \\
	-O b \\
	-o TMP/${prefix}.bcf \\
	'${vcf}'

mv  TMP/${prefix}.bcf ./
mv  TMP/${prefix}.bcf.csi ./

cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
}
