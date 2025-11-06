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
process SAMTOOLS_SAMPLES {
tag "${meta.id?:""}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path("REFS/*")
    tuple val(meta ),path("BAMS/*")
output:
    tuple val(meta),path("*.tsv"),emit:samplesheet
    path("versions.yml"),emit:versions
script:
    def args1 = task.ext.args1?:""
    // check sample is defined
    def check_sample = (task.ext.check_sample?:true).toBoolean()
    // check no dup sample
    def check_dup = (task.ext.check_dup?:true).toBoolean()
    // check fasta is known
    def check_ref = (task.ext.check_ref?:true).toBoolean()
"""
mkdir -p TMP

# not created if there is no ref...
mkdir -p REFS

find \${PWD}/BAMS/ -name "*am" |\\
    samtools samples \\
        ` find REFS/ \\( -name "*.fasta" -o  -name "*.fa" -o -name "*.fna" \\)  -printf "-f %p "` \\
        ${args1} > TMP/samplesheet.tsv

if ${check_sample}
then
	awk -F '\t' '(\$1==".")' TMP/samplesheet.tsv > TMP/sn.txt
	sed 's/^/MISSING RG:SM /' TMP/sn.txt 1>&2
	test ! -s TMP/sn.txt
fi

if ${check_ref}
then
	awk -F '\t' '(\$3==".")' TMP/samplesheet.tsv > TMP/fa.txt
	sed 's/^/MISSING REF /' TMP/fa.txt 1>&2
	test ! -s TMP/fa.txt
fi

if ${check_dup}
then

	cut -f1 TMP/samplesheet.tsv | sort | uniq -d > TMP/dups.txt
	sed 's/^/DUPLICATE SAMPLE: /' TMP/dups.txt 1>&2
	test ! -s TMP/dups.txt

	cut -f2 TMP/samplesheet.tsv | sort | uniq -d > TMP/dups.txt
	sed 's/^/DUPLICATE BAM: /' TMP/dups.txt 1>&2
	test ! -s TMP/dups.txt

fi



mv TMP/samplesheet.tsv ./



cat << END_VERSIONS > versions.yml
"${task.process}":
	samtools: "\$(samtools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""

stub:
"""
find \${PWD}/BAMS/ -name "*am" | sort -V | awk '{printf("S%d\t%s\t%s\\n",NR,\$0,"${params.fasta?:"."}");}'  > samplesheet.tsv
touch versions.yml
"""
}
