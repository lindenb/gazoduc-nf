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
process FASTA_TO_CONTIGS_01 {
tag "${reference}"
executor "local"
input:
	val(meta)
	val(reference)
output:
	path("contigs.tsv"),emit:output
	path("contigs.bed"),emit:bed
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
test -f "${reference}.fai"

awk -F '\t' 'BEGIN{printf("contig\tsize\tindex\treference\\n");} {printf("%s\t%s\t%d\t${reference}\\n",\$1,\$2,NR);}' '${reference}.fai' > contigs.tsv
awk -F '\t' '{printf("%s\t0\t%s\t${reference}\\n",\$1,\$2);}' '${reference}.fai' | LC_ALL=C sort -t '\t' -k1,1 -k2,2n > contigs.bed


##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">get contigs from reference</entry>
        <entry key="reference">${reference}</entry>
</properties>
EOF
"""
}
