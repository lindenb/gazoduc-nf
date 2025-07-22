

process APPLY_INDEXCOV {
tag "${meta.id}"
afterScript "rm -rf TMP"
label "process_high"
conda "${moduleDir}/../../../conda/goleft.yml"
when:
        task.ext.when == null || task.ext.when
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta),path("BAMS/*")
output:
    tuple val(meta),path ("*indexcov.zip"),emit: zip
    tuple val(meta),tuple path("*.indexcov.bed.gz"), path("*.indexcov.bed.gz.tbi"),emit:bed
    path "versions.yml" , emit: versions
script:
    def name = meta.id
    def prefix = task.ext.prefix?:name
    def prefix2 = prefix +"."
"""
hostname 1>&2
mkdir -p "${prefix}" TMP

find BAMS/  -name "*.bam" -o -name "*.cram" > TMP/jeter.list

# test not empty
test -s TMP/jeter.list

# test all files exist
xargs -a TMP/jeter.list -L1 --verbose test -f

sed 's/\\.cram/.cram.crai/' TMP/jeter.list > TMP/tmp.bams.list

# test all files exist
xargs -a TMP/tmp.bams.list -L1 --verbose test -f


goleft indexcov \\
    --fai "${fasta}.fai"  \\
    --excludepatt `awk -F '\t' '(!(\$1 ~ /^(chr)?[0-9XY]+\$/)) {print \$1;}' "${fai}" | paste -s -d '|' `  \\
    --sex `awk -F '\t' '(\$1 ~ /^(chr)?[XY]\$/) {print \$1;}' "${fai}" | paste -s -d, ` \\
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
}
