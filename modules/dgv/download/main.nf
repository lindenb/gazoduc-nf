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

/**
DGV (Database of Genomic Variants) as BED file from ${url}
*/
process DGV_DOWNLOAD {
label "process_single"
tag "${meta.id}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta),path(dict)
output:
    tuple val(meta),path("*.bed.gz"),path("*.bed.gz.tbi"),emit:bed
    path("versions.yml"),emit:versions
script:
	def url= task.ext.url?:""
	def prefix = task.ext.prefix?:"${dict.baseName}.dgv"
	def awk_expr = task.ext.awk_expr?:"(1==1)"
	def jvm = ""

	if("${url}"=="") {
		if(meta.ucsc_name=="hg38") {
			url = "https://dgv.tcag.ca/dgv/docs/GRCh38_hg38_variants_2025-12-01.txt"
			}
		else if(meta.ucsc_name=="hg19") {
			url = "https://dgv.tcag.ca/dgv/docs/GRCh37_hg19_variants_2025-12-01.txt"
			}
		else
			{
			log.warn("${task.proces} : undefined ucsc_name")
			}
		}
"""
mkdir TMP

if test ! -z "${url}"
then

	curl -s -L "${url}" |\\
		awk -F '\t' '(${awk_expr}){printf("%s\t%s\t%s\t%s",\$2,\$3,\$4,\$1);for(i=5;i<=NF;i++){printf("\t%s",\$i);}printf("\\n");}' |\\
		sed  's/^chr/#chr/' |\\
		jvarkit ${jvm} bedrenamechr -f "${dict}" --column 1 --convert SKIP |\\
		LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n > TMP/dgv.bed

else
	touch TMP/dgv.bed
fi


bgzip --force TMP/dgv.bed
tabix --comment '#'  -p bed -f TMP/dgv.bed.gz

mv TMP/dgv.bed.gz ./${prefix}.bed.gz
mv TMP/dgv.bed.gz.tbi ./${prefix}.bed.gz.tbi

cat << EOF > versions.yml
${task.process}:
	url: ${url}
EOF
"""


stub:
def prefix = task.ext.prefix?:"${dict.baseName}.dgv"
"""
touch versions.yml ${prefix}.bed.gz ${prefix}.bed.gz.tbi
"""
}
