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
process HC_GENOMICDB_IMPORT {
tag "${bed.name}"
label "process_single"
afterScript "rm -rf TMP"
memory = {20.GB  * task.attempt}
errorStrategy "retry"
maxRetries 2
time "3h"
input:
        tuple val(meta1),path(fasta)
        tuple val(meta2),path(fai)
        tuple val(meta3),path(dict)
        tuple path(bed),path("VCFS/*")
output:
        tuple path(bed),path("database"),emit:output
script:
	def batchSize=-1
"""
hostname 1>&2
module load gatk/0.0.0
mkdir -p TMP
set -x

find ./VCFS -name "*.vcf.gz" > TMP/jeter.list
test -s TMP/jeter.list

cat TMP/jeter.list | while read F
do
        gunzip -c "\${F}" | grep "^#CHROM" -m1 | cut -f 10 | tr "\\n" "\t" >> TMP/sample.map
        echo "\${F}" >> TMP/sample.map
done


# sort on sample name
LC_ALL=C sort -t '\t' -k1,1 -T TMP TMP/sample.map > "TMP/jeter.map"
test -s TMP/jeter.map
mv TMP/jeter.map TMP/sample.map



SQRT=`awk 'END{X=NR;if(X<10){print(X);} else {z=sqrt(X); print (z==int(z)?z:int(z)+1);}}' TMP/sample.map`

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" GenomicsDBImport \\
    -R ${fasta} \\
    --batch-size ${(batchSize as Integer) <= 0 ? "\${SQRT}" : ""+batchSize} \\
    ${(task.cpus as Integer) > 1 ? "  --reader-threads " +task.cpus : "" } \\
    --sample-name-map TMP/sample.map \
    -L "${bed}" \
    --genomicsdb-workspace-path "TMP/database"


mv TMP/database ./
"""
}
