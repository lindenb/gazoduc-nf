
process BEDTOOLS_MAKEWINDOWS {
tag "${bed.name}"
label "process_quick"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta),path(bed)
output:
    tuple val(meta),path("*.bed"),emit:bed
    path("versions.yml"),emit:versions
script:
    def option1= (bed.name.endsWith(".fai")?"-g":"-b")
    def args = task.ext.args?:""
    if(args.trim().isEmpty()) throw new IllegalArgumentException("args empty for ${task.process}")
    def prefix = task.ext.prefix?:bed.baseName+".makewindows"
"""
mkdir -p TMP
bedtools makewindows ${option1} "${bed}" ${args} |\\
    sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter.bed

mv TMP/jeter.bed "${prefix}.bed"


cat << END_VERSIONS > versions.yml
"${task.process}":
    bcftools: \$(bedtools --version | awk '(NR==1)  {print \$NF}')
END_VERSIONS
"""
}
