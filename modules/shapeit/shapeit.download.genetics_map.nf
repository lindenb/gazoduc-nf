/*

Copyright (c) 2023 Pierre Lindenbaum

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


process SHAPEIT_DOWNLOAD_GENETICS_MAP {
tag "${genomeId}"
afterScript "rm -f jeter.tsv jeter.tar.gz"
input:
	val(meta)
	val(genomeId)
output:
	path("${genomeId}.genetics_map.tsv"),emit:output // file with header and columns contig/path/genome_id
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def url = genome.shapeit_genetic_map
"""
hostname 1>&2
test ! -z "${url}"

wget -O jeter.tar.gz "${url}"
tar xvfz jeter.tar.gz
rm jeter.tar.gz

find \${PWD}/ -type f -name "*.gmap.gz" |\
	awk -F '/' 'BEGIN{printf("contig\tgmap\tgenomeId\\n");} {F=\$NF;split(F,a,/[\\.]/); printf("%s\t%s\t${genomeId}\\n",a[1],\$0);}' > jeter.tsv

mv jeter.tsv "${genomeId}.genetics_map.tsv"

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">download genetics map for shapeit</entry>
        <entry key="url"><a>${url}</a></entry>
</properties>
EOF
"""
}
