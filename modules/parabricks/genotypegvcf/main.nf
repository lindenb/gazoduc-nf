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
process PB_GENOTYPEGVCF {
  tag "${meta.id?:""}"
  label 'process_gpu'
  afterScript "rm -rf TMP"
  input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta),path(gvcf),path(gvcf_idx)
  output:
	tuple val(meta),
		path("${meta.id}.vcf.gz"),
		path("${meta.id}.vcf.gz.tbi"),emit:bams
	path("${meta.id}.hc.log"),emit:log
	path("versions.yml"),emit:versions
  script:
    def prefix = task.ext.prefix?:"${meta.id}"
 """
mkdir -p TMP

nvidia-smi 1>&2

pbrun genotypegvcf \\
    --num-gpus ${task.ext.gpus} \\
    --ref ${fasta} \\
    --in-gvcf "${gvcf}" \\
    --out-vcf "${prefix}.vcf.gz" \\
    --tmp-dir TMP \\
    --logfile ${prefix}.log


cat <<-END_VERSIONS > versions.yml
"${task.process}":
        pbrun: \$(grep Parabr -m1 -A1 "${prefix}.log" | grep -o 'Version [^ ]*' )
END_VERSIONS
"""

 stub:
    def sample = meta.id
    """
    touch ${sample}.g.vcf.gz
    touch ${sample}.g.vcf.gz.tbi
    touch ${sample}.hc.log
    touch versions.yml
    """
}
