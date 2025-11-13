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
process PB_DEEPVARIANT {
  tag "${meta.id}"
  label 'process_gpu'
  label 'parabricks'

  afterScript "rm -rf TMP"
  input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta),path(bam),path(bai)
  output:
	  tuple val(meta), path("${meta.id}.g.vcf.gz"),path("${meta.id}.g.vcf.gz.tbi"),emit:gvcf
    path("${meta.id}.deepvariant.log")
    path("versions.yml"),emit:versions
  script:
    def sample = meta.id
    def args1 = task.ext.args1?:""
    def num_gpus = task.accelerator ? "--num-gpus $task.accelerator.request" : ''
 """
	mkdir -p TMP

  nvidia-smi 1>&2
  set -x
  pbrun deepvariant \\
      ${num_gpus} \\
      --ref ${fasta} \\
      --in-bam "${bam}" \\
      --gvcf \\
      --out-variants "${sample}.g.vcf.gz" \\
      --tmp-dir TMP \\
      --logfile ${sample}.deepvariant.log \\
      ${args1}


cat << END_VERSIONS > versions.yml
"${task.process}":
    pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
END_VERSIONS
"""


stub:
def sample = meta.id
"""
  touch ${sample}.g.vcf.gz
  touch ${sample}.g.vcf.gz.tbi
  touch ${sample}.log
  touch ${sample}.deepvariant.log
  touch versions.yml
"""
}
