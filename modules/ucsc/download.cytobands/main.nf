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
include {isBlank} from '../../utils/functions.nf'

process DOWNLOAD_CYTOBAND {
tag "${meta.id}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta),path(dict)
output:
	tuple val(meta),path("*.bed.gz"),path("*.bed.gz"),emit:tabix
	tuple val(meta),path("*.bed"),emit:bed
	path("versions.yml"),emit:versions
script:	
	def ucsc_name = (meta.ucsc_name?:"")
	def url = task.ext.url?:""
	if(task.ext.url==null) {
		url = "https://hgdownload.cse.ucsc.edu/goldenPath/${ucsc_name}/database/cytoBandIdeo.txt.gz"	
		}
	def prefix = task.ext.prefix?:"${meta.id}"

"""
hostname 1>&2
mkdir -p TMP


if ${isBlank(url)}
then

	touch "${prefix}.cytoBandIdeo.bed"
	

else

	curl -L  "${url}" |\\
			gunzip -c |\\
			jvarkit bedrenamechr -f "${dict}" --column 1 --convert SKIP |\\
			LC_ALL=C sort -S ${task.memory.kilo}  -T TMP -t '\t' -k1,1 -k2,2n  > "${prefix}.cytoBandIdeo.bed" 

fi


cat "${prefix}.cytoBandIdeo.bed" | bgzip >  "${prefix}.cytoBandIdeo.bed.gz"
tabix -p bed   "${prefix}.cytoBandIdeo.bed.gz"

cat << EOF > versions.yml
"${task.process}":
	url: "${url}"
EOF
"""

stub:
	def prefix = task.ext.prefix?:"${meta.id}"
"""
touch versions.yml "${prefix}.cytoBandIdeo.bed" "${prefix}.cytoBandIdeo.bed.gz" "${prefix}.cytoBandIdeo.bed.gz.tbi"
"""
}

