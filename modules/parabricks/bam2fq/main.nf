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
process PB_BAM2FQ {
  tag "${meta.id?:""}"
  label 'process_gpu'
  label 'parabricks'
  afterScript "rm -rf TMP"
  input:
    tuple val(meta),path(bam),path(fasta),path(fai)
  output:
      tuple val(meta), path("*.fastq.gz", arity: '1..*'), emit: fastqs
      path("${meta.id}.bam2fq.log")
      path("versions.yml"),emit:versions
  script:
    def sample = meta.id?:bam.simpleName
    def with_orphan = task.ext.with_orphan?:false
    def with_singleton = task.ext.with_singleton?:false
 """
	mkdir -p TMP/OUT

  pbrun bam2fq \\
      --num-threads ${task.cpus} \\
      --ref ${fasta} \\
      --out-prefix "TMP/OUT/${sample}." \\
      --out-suffixF  R1.fastq.gz \\
      --out-suffixF2 R2.fastq.gz \\
      ${with_orphan   ?"--out-suffixO  .orphan1.fastq.gz":""} \\
      ${with_orphan   ?"--out-suffixO2 .orphan2.fastq.gz":""} \\
      ${with_singleton?"--out-suffixS  .singletons.fastq.gz":""} \\
      --in-bam "${bam}" \\
      --tmp-dir TMP \\
      --logfile ${sample}.bam2fq.log

find .  1>&2

mv -v TMP/OUT/*.gz ./

cat << END_VERSIONS > versions.yml
"${task.process}":
    pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
END_VERSIONS
"""


stub:
    def sample = meta.id
    """
    touch ${sample}.R1.fastq.gz
    touch ${sample}.R2.fastq.gz
    touch ${sample}.bam2fq.txt
    touch versions.yml
    """
}
