/*

Copyright (c) 2024 Pierre Lindenbaum

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
process HC_COMBINE2 {
tag "${bed.name}"
label "process_quick"
afterScript "rm -rf TMP"
input:
        tuple val(meta1),path(fasta)
        tuple val(meta2),path(fai)
        tuple val(meta3),path(dict)
	tuple val(meta4),path(dbsnp)
	tuple val(meta5),path(dbsnp_tbi)
        tuple path(bed),path("VCFS/*")
output:
        tuple path(bed),path("*.g.vcf.gz"),path("*.g.vcf.gz.tbi"),emit:output
script:
"""
hostname 1>&2
module load gatk/0.0.0
mkdir -p TMP
set -x
MD5=`cat TMP/jeter.list | sort | sha1sum | cut -d ' ' -f1`


find VCFS/ -name "*.g.vcf.gz" | sort > TMP/jeter.list

        gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" \\
                CombineGVCFs \\
                -R "${fasta}" \\
                -L "${bed}" \\
                -V TMP/jeter.list \\
                -O "TMP/combine2.g.vcf.gz" \\
                -G StandardAnnotation \\
                -G AS_StandardAnnotation


mv TMP/combine2.g.vcf.gz ./\${MD5}.g.vcf.gz
mv TMP/combine2.g.vcf.gz.tbi ./\${MD5}.g.vcf.gz.tbi
"""
}
