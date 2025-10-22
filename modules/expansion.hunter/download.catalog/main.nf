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
process XHUNTER_DOWNLOAD_CATALOG {
label "process_single"
tag "${meta.ucsc_name}"
input:
	tuple val(meta ),path(fasta)
	tuple val(meta2),path(fai)
output:
	tuple val(meta),path("*.json"),optional:true,emit:catalog
	path("versions.yml"),emit:versions
script:
	def url0  = (meta.ucsc_name==null || !meta.ucsc_name.matches("hg(19|38)"))?"":"https://github.com/Illumina/RepeatCatalogs/blob/master/${meta.ucsc_name}/variant_catalog.json?raw=true"
	def url = task.ext.url?:url0
if(url.isEmpty())
"""
touch versions.yml
"""
else
"""
curl -L -o variant_catalog.${meta.ucsc_name}.json "${url}"

cat << EOF > versions.yml
"${task.process}":
	url : "${url}"
EOF
"""

stub:
"""
touch versions.yml ${meta.id}.catalog.json
"""
}
