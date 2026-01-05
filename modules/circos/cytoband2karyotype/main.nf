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
FOR GLICI users


note to self: sed -i 's/p@univ/p\\@univ/g' (..)/lib/perl5/5.32/core_perl/Config_heavy.pl 

*/

process CYTOBAND_TO_KARYOTYPE {
tag "${meta.id}"
label "process_single"
conda "${moduleDir}/../../../conda/circos.yml"
afterScript "rm -rf TMP"
input:
    tuple val(meta),path(fai)
    tuple val(meta2),path(cytoband)//optional
output:
	tuple val(meta),path("*.txt"),emit:karyotype
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix ?: "${meta.id}.karyotype"
    def args1 = task.ext.args1?:""
"""
mkdir -p TMP

awk '{
    C=(NR%2==0?"200,200,200":"100,100,100");
    if(C ~ /^(chr)?[0-9XY]+\$/) C=\$1;
    printf("chr - %s %s 0 %s %s\\n",\$1,\$1,\$2, C);}
    ' "${fai}" > TMP/jeter.txt

# chrX    42400000        46400000        p11.3   gpos75
# chrX    46400000        49800000        p11.23  gneg
# chrX    49800000        54800000        p11.22  gpos25
# chrX    54800000        58100000        p11.21  gneg

if ${cytoband?true:false}
then
    awk -F '\t' '(\$4!="" && \$5!="") {
        printf("band %s %s %s %s %s %s\\n",\$1,\$4,\$4,\$2,\$3,\$5);
        }' "${cytoband}" >> TMP/jeter.txt
fi

mv TMP/jeter.txt ${prefix}.txt

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    awk: todo
END_VERSIONS
"""

stub:
"""
touch versions.yml "${meta.id}.txt"
"""
}
