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
/**
ORIGINAL snakemake workflow by Raphael Blanchet PhD.


*/
include { parseBoolean} from '../../../modules/utils/functions.nf'

process PB_FQ2BAM {
  tag "${meta.id}"
  label 'process_gpu'
  label 'parabricks'
  afterScript "rm -rf TMP"

  input:
    tuple val(meta2),path(fasta)
    tuple val(meta3),path(fai)
    tuple val(meta4),path(bwa_index_dir)
    tuple val(meta5),path(known_indels), path(known_indels_tbi)
    tuple val(meta), path(fq_R1), path(fq_R2)
  output:
	    tuple val(meta), path("*.cram") ,path("*.crai"),emit: bam,    optional: true
      tuple val(meta), path("*.duplicate.metrics.txt"), emit:duplicate_metrics, optional:true
      tuple val(meta), path("*.bqsr.report.txt"), emit:bqsr_report, optional:true
      tuple val(meta), path("*.qc_metrics"),optional:true, emit:qc_metrics
      tuple val(meta), path("*.md5"),optional:true, emit:md5
      path("*.log")
      path("versions.yml"),emit:versions
  script:

    def sample = task.ext.prefix?:meta.id
    def lib = meta.LIB?:sample
    def rgid = meta.rgid?:sample
    def pl = meta.PL?:"ILLUMINA"
    def with_bqsr = parseBoolean(task.ext.with_bqsr?:false)
    def cpu_per_gpu = task.ext.cpu_per_gpu?:1
    def bqsr_args = with_bqsr && known_indels? "--knownSites \"${known_indels.name}\"   --out-recal-file \"${sample}.bqsr.report.txt\"   " : ""
    def low_memory = task.ext.low_memory==true || task.attempt>1 ? "--low-memory" : ""
    def fixmate_args = (parseBoolean(task.ext.with_fixmate?:true)?"--fix-mate":"")
    //TODO   add    ‑‑out‑qc‑metrics‑dir "TMP/OUT/${sample}.qc_metrics" 
 """
	mkdir -p TMP/TMP
	mkdir -p TMP/REF
  mkdir -p TMP/OUT
  
	# pb doesn't like the symlinks ?

  find ${bwa_index_dir}/  1>&2
  cp -v  ${bwa_index_dir}/*.ann TMP/REF/${fasta.name}.ann
  cp -v  ${bwa_index_dir}/*.pac TMP/REF/${fasta.name}.pac
  cp -v  ${bwa_index_dir}/*.sa  TMP/REF/${fasta.name}.sa
  cp -v  ${bwa_index_dir}/*.amb TMP/REF/${fasta.name}.amb
  cp -v  ${bwa_index_dir}/*.alt TMP/REF/${fasta.name}.alt || true
  cp -v  ${bwa_index_dir}/*.bwt TMP/REF/${fasta.name}.bwt
  cp -v  "${fasta}"  TMP/REF/${fasta.name}
  cp -v  "${fai}"  TMP/REF/${fasta.name}.fai
	ls -lah 1>&2

	# show what's here
	find .  1>&2
	pwd 1>&2

  nvidia-smi 1>&2
   
  pbrun fq2bam \\
      --num-gpus ${task.ext.gpus} \\
      --gpusort \\
      --bwa-cpu-thread-pool ${cpu_per_gpu} \\
      --num-cpu-threads-per-stage ${cpu_per_gpu} \\
      --memory-limit ${task.memory.giga} \\
      --ref TMP/REF/${fasta.name} \\
      --in-fq ${fq_R1} ${fq_R2?"${fq_R2}":""} \\
      --out-bam "TMP/OUT/${sample}.cram" \\
      --read-group-sm "${sample}" \\
      --read-group-lb "${lib}" \\
      --read-group-pl "${pl}" \\
      --read-group-id-prefix "${rgid}" \\
      --gpusort \\
      --out-duplicate-metrics "TMP/OUT/${sample}.duplicate.metrics.txt" \\
      ${bqsr_args} \\
      ${fixmate_args} \\
      --tmp-dir TMP \\
      --logfile ${sample}.log \\
      ${low_memory} 1>&2

find .  1>&2

md5sum < "TMP/OUT/${sample}.cram" > "TMP/OUT/${sample}.cram.md5"

mv "TMP/OUT/${sample}.duplicate.metrics.txt" ./ || true
mv "TMP/OUT/${sample}.qc_metrics" ./ || true
mv "TMP/OUT/${sample}.cram" ./
mv "TMP/OUT/${sample}.cram.crai" ./
mv "TMP/OUT/${sample}.cram.md5" ./

cat << END_VERSIONS > versions.yml
"${task.process}":
    pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
END_VERSIONS
"""


stub:
    def sample = meta.id
    """
    touch ${sample}.cram
    touch ${sample}.cram.md5
    touch ${sample}.cram.crai
    touch ${sample}.duplicate.metrics.txt
    touch versions.yml
    """
}
