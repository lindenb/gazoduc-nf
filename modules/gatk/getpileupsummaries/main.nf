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

process GET_PILEUP_SUMMARIES {
tag "${meta.id?:""}"
label "process_single"
afterScript 'rm -rf TMP'
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta4 ),path(vcf),path(vcfidx)
    tuple val(meta5 ),path(opt_file),path(opt_file_idx)//whatever kind of file used for option '-L'
	tuple val(meta ),path(bam),path(bai)
output:
    tuple val(meta),path(".table"),emit:table
    path("versions.yml"),emit:versions
when:
    task.ext.when == null || task.ext.when
script:
    def jvm = task.ext.jvm?:"-Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP"
"""
hostname 1>&2
mkdir -p TMP

gatk --java-options "${jvm}" GetPileupSummaries \\
    -I ${bam} \\
    -V ${vcf} \\
    -L ${opt_file?opt_file:vcf} \\
    -O TMP/jeter.pileups.table

mv TMP/jeter.pileups.table "${meta.id}.pileups.table"

cat << EOF > version.yml
${task.process}:
    gatk: "\$( gatk --version 2> /dev/null  | paste -s -d ' ' )"
EOF
"""
stub:
"""
touch versions.yml "${meta.id}.pileups.table"
"""
}
