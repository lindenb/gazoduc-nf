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

process APPLY_INDEXCOV {
tag "${meta.id}"
afterScript "rm -rf TMP"
label "process_high"
conda "${moduleDir}/../../../conda/goleft.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta),path("BAMS/*")
output:
    tuple val(meta),path ("*indexcov.zip"),emit: zip
    tuple val(meta),path("*.indexcov.bed.gz"), path("*.indexcov.bed.gz.tbi"),emit:bed
    path "versions.yml" , emit: versions
script:
    def prefix = task.ext.prefix?:"${meta.id}"
    def prefix2 = prefix +"."
"""
hostname 1>&2
mkdir -p "${prefix}" TMP

find BAMS/ \\( -name "*.bam" -o -name "*.cram" \\) | sort -V -T TMP > TMP/jeter.list

# test not empty
test -s TMP/jeter.list

# test all files exist
xargs -a TMP/jeter.list -L1 --verbose test -f

sed 's/\\.cram/.cram.crai/' TMP/jeter.list | sort | uniq > TMP/tmp.bams.list

# test all files exist
xargs -a TMP/tmp.bams.list -L1 --verbose test -f

goleft indexcov \\
    --fai "${fasta}.fai"  \\
    `awk -F '\t' '(!(\$1 ~ /^(chr)?[0-9XY]+\$/)) {print \$1;}' "${fai}" | paste -s -d '|'  | awk '(length(\$0)>0)  {printf("--excludepatt %s",\$0);}'`  \\
    `awk -F '\t' '(\$1 ~ /^(chr)?[XY]\$/) {print \$1;}' "${fai}" | paste -s -d, | awk '(length(\$0)>0)  {printf("--sex %s",\$0);}' ` \\
    --directory "${prefix}" \\
    `awk '/.crai\$/ {X=1;} END {if(X==1) printf(" --extranormalize ");}' TMP/tmp.bams.list` \\
    `cat TMP/tmp.bams.list`


#create tabix index
tabix -f -p bed "${prefix}/${prefix}-indexcov.bed.gz"

ln -s "${prefix}/${prefix}-indexcov.bed.gz" "${prefix}.indexcov.bed.gz"
ln -s "${prefix}/${prefix}-indexcov.bed.gz.tbi" "${prefix}.indexcov.bed.gz.tbi"



# all the generated files
find "${prefix}" -type f > "${prefix2}files.list"

# zip results
cat "${prefix2}files.list" | zip -@ -9  "${prefix2}indexcov.zip"


cat << EOF > versions.yml
"${task.process}":
    indexcov: todo
EOF
"""

stub:
    def prefix = task.ext.prefix?:"${meta.id}"
    def prefix2 = prefix +"."
"""
mkdir -p indexcov
touch "${prefix}.indexcov.bed.gz"
touch "${prefix}.indexcov.bed.gz.tbi"
touch indexcov/"${prefix}-indexcov.bed.gz"
touch indexcov/"${prefix}-indexcov.bed.gz.tbi"


zip -r "${prefix2}indexcov.zip" indexcov
touch versions.yml 
"""
}
