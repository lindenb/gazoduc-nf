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
include {k1_signature} from '../../../modules/utils/k1.nf'


workflow ALPHAMISSENSE {
	take:
		meta
		fasta
		fai
		dict
		bed
		vcfs /* meta, vcf,vcf_index */
	main:
		versions = Channel.empty()
		
		DOWNLOAD(fasta,fai,dict,bed)
		versions = versions.mix(DOWNLOAD.out.versions)

		DOWNLOAD.out.bed.view{"it1 ${it}"}
		DOWNLOAD.out.tbi.view{"it2 ${it}"}
		DOWNLOAD.out.header.view{"it3 ${it}"}
		vcfs.view{"it4 ${it}"}

		ANNOTATE(DOWNLOAD.out.bed, DOWNLOAD.out.tbi, DOWNLOAD.out.header, vcfs)
		versions = versions.mix(ANNOTATE.out.versions)
	emit:
		vcf = ANNOTATE.out.vcf
		versions
	}


process DOWNLOAD {
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
	tuple val(meta4),path(bed)
output:
	tuple val(meta1),path("*.tsv.gz"),emit:bed
	tuple val(meta1),path("*.tsv.gz.tbi"),emit:tbi
	tuple val(meta1),path("*.header"),emit:header
	tuple val(meta1),path("*.md"),emit:doc
	path("versions.yml"),emit:versions
script:
    def k1 = k1_signature()
    def TAG = task.ext.tag?:"ALPHAMISSENSE"
    def whatis="Data from AlphaMissense https://www.science.org/doi/10.1126/science.adg7492 Accurate proteome-wide missense variant effect prediction with AlphaMissense"
    def base="https://storage.googleapis.com/dm_alphamissense/AlphaMissense_"
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP

cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg38}\t${base}hg38.tsv.gz
1:${k1.hg19}\t${base}hg19.tsv.gz
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' |\\
	sed 's/^chr//' |\\
	sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv

join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv |\\
	sort |\\
	uniq > TMP/jeter.url

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


cat << EOF > ${TAG}.md
${whatis}
EOF

cat << EOF > versions.yml
${task.process}:
	jvarkit: TODO
EOF
"""
}


process ANNOTATE {
tag "${meta.id?:vcf.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(tabix)
    tuple val(meta2),path(tbi)
    tuple val(meta3),path(header)
    tuple val(meta ),path(vcf),path(vcf_idx)
output:
    tuple val(meta),path("*.bcf"),path("*.csi"),emit:vcf
	path("versions.yml"),emit:versions
script:
   	def TAG = task.ext.tag?:"ALPHAMISSENSE"
	def prefix = task.ext.prefix?:vcf.baseName+".alphamissense"
"""
hostname 1>&2
mkdir -p TMP OUTPUT

bcftools annotate \\
    --threads ${task.cpus} \\
    -a "${tabix}" \\
    -h "${header}" \\
    -c "CHROM,POS,REF,ALT,${TAG}_PATHOGENOCITY,${TAG}_CLASS" \\
    --merge-logic "${TAG}_PATHOGENOCITY:max,${TAG}_CLASS:unique" \\
    -O b -o TMP/${prefix}.bcf '${vcf}'
 
bcftools index \\
    --threads ${task.cpus} \\
    -f TMP/${prefix}.bcf

mv TMP/*.bcf ./
mv TMP/*.bcf.csi ./

cat << EOF > versions.yml
${task.process}:
	jvarkit: TODO
EOF
"""
}
