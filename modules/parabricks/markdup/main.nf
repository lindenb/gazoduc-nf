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
process PB_MARKDUP {
  tag "${meta.id?:""}"
  label 'process_gpu'
  label 'parabricks'
  afterScript "rm -rf TMP"
  input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta),path(bam)
  output:
	tuple val(meta),
		path("*.bam"),
		path("*.bai"),emit:bam
    tuple val(meta),path("${prefix}.metrics"),emit:metrics
	path("${meta.id}.hc.log"),emit:log
	path("versions.yml"),emit:versions
  script:
    def args1  = task.args1?:"" //e.g: --markdups-assume-sortorder-queryname
    def prefix = task.ext.prefix?:"${meta.id}.markdup"
 """
mkdir -p TMP

nvidia-smi 1>&2

pbrun markdup \\
    --num-gpus ${task.ext.gpus} \\
    --ref ${fasta} \\
    --in-bam "${bam}" \\
    --out-bam  "TMP/${prefix}.bam" \\
    ‑‑num‑worker‑threads ${task.cpus} \\
    ‑‑out‑duplicate‑metrics ${prefix}.metrics \\
    --tmp-dir TMP \\
    --logfile ${prefix}.hc.log \\

mv "TMP/${prefix}.bam" ./
mv "TMP/${prefix}.bam.bai" ./

cat <<-END_VERSIONS > versions.yml
"${task.process}":
        pbrun: \$(grep Parabr -m1 -A1 "${sample}.hc.log" | grep -o 'Version [^ ]*' )
END_VERSIONS
"""

 stub:
    def sample = meta.id
    """
    touch ${sample}.bam
    touch ${sample}.baam.bai
    touch ${sample}.hc.log
    touch versions.yml
    """
}
