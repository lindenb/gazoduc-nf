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

process DOWNLOAD_REFGENE {
tag "${meta.id}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta),path(dict)
output:
	tuple val(meta),path("*.txt.gz"),path("*.txt.gz.tbi"),emit:tabix
	path("versions.yml"),emit:versions
script:
	def ucsc_name = (meta.ucsc_name?:"")
	def url = task.ext.url?:""
	if(task.ext.url==null) {
		url = "https://hgdownload.cse.ucsc.edu/goldenPath/${ucsc_name}/database/refGene.txt.gz"	
		}
	def prefix = task.ext.prefix?:"${meta.id}"
"""
hostname 1>&2
mkdir -p TMP
if ${!url.isEmpty()}
then

curl -L "${url}" |\\
		gunzip -c |\\
		jvarkit bedrenamechr -f "${dict}" --column 3 --convert SKIP |\\
		LC_ALL=C sort -S ${task.memory.kilo}  -T TMP -t '\t' -k3,3 -k5,5n |\\
		bgzip > "TMP/${prefix}.refGene.txt.gz"

else

	touch TMP/${prefix}.refGene.txt
	bgzip TMP/${prefix}.refGene.txt

fi

tabix -f -0 -b 5 -e 6 -s 3 "TMP/${prefix}.refGene.txt.gz"

mv "TMP/${prefix}.refGene.txt.gz" ./
mv "TMP/${prefix}.refGene.txt.gz.tbi" ./

cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${url}"
END_VERSIONS
"""

stub:
	def prefix = task.ext.prefix?:"${meta.id}"
"""
touch versions.yml "${prefix}.refGene.txt.gz" "${prefix}.refGene.txt.gz.tbi" 
"""
}

