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

/**
 *
 * builds a file to convert the chrX to X with bcftools annotate for
 * rvtest because the software is buggy...
 * https://github.com/zhanxw/rvtests/issues/80
 *
 **/
process RVTESTS_REHEADER_01 {
tag "${reference}"
executor "local"
input:
	val(meta)
	val(reference)
output:
	path("reheader.tsv"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2

cut -f 1 "${reference}.fai" > chroms.A.txt
sed 's/^chr//' chroms.A.txt > chroms.B.txt
paste chroms.A.txt chroms.B.txt > reheader.tsv
test -s reheader.tsv
rm chroms.A.txt chroms.B.txt

############################################
cat << EOF > version.xml
<properties id="${task.process}">
  <entry key="name">${task.process}</entry>
  <entry key="description">creates a chrom notation converter for bcftools annotate / rvtests.  <url>https://github.com/zhanxw/rvtests/issues/80</url>.</entry>
</properties>
EOF
"""
}
