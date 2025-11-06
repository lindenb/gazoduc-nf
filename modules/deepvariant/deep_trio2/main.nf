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

process DEEP_TRIO2 /* 1 child+1 parent */{
label "process_single"
tag "${meta.id} ${C_bam.name} ${optional_bed?optional_bed.name:""}"
afterScript "rm -rf TMP  .keras  .parallel"
conda "${moduleDir}/../../../conda/deepvariant.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta ),path(C_bam),path(C_bai), // CHILD
	 	path(P_bam),path(P_bai), //parent
		path(optional_bed)
output:
	tuple val(meta),
		path("*.C.g.vcf.gz"),path("*.C.g.vcf.gz.tbi"),
		path("*.P.g.vcf.gz"),path("*.P.g.vcf.gz.tbi"),
		path(optional_bed),optional:true,emit:gvcf
	
	tuple val(meta),
		path("*.C.std.vcf.gz"),path("*.C.std.vcf.gz.tbi"),
		path("*.P.std.vcf.gz"),path("*.P.std.vcf.gz.tbi"),
		path(optional_bed),optional:true,emit:vcf
	
	path("versions.yml"),emit:versions
script:
	if(!meta.parent) throw new IllegalArgumentException("${task.process} missing meta.parent")
	def model_type = task.ext.model_type?:(meta.model_type?:"")
	if(model_type.isEmpty()) throw new IllegalArgumentException("${task.process} missing ext.model_type (e.g WGS,WES,PACBIO,ONT_R104,HYBRID_PACBIO_ILLUMINA)");
	def child = meta.id
	def parent = meta.parent
	def prefix = task.ext.prefix?:"deepvariant"
	def keep_gvcf= (task.ext.keep_gvcf?:true).toBoolean()
	def keep_vcf= (task.ext.keep_vcf?:true).toBoolean()
"""
	hostname 1>&2
	mkdir -p TMP/TMP
	
	run_deeptrio \\
		--model_type ${model_type} \\
		--ref "${fasta}" \\
		--reads_child "${C_bam}" \\
		--reads_parent "${P_bam}" \\
		--sample_name_child "${child}" \\
	    --sample_name_parent1 "${father}" \\
		--num_shards ${task.cpus} \\
		--intermediate_results_dir TMP/TMP \\
		${optional_bed?"--regions ${optional_bed}":""} \\
		--output_vcf_child TMP/child.vcf.gz \\
  		--output_vcf_parent1 TMP/parent.vcf.gz \\
  		--output_gvcf_child  TMP/child.g.vcf.gz \\
  		--output_gvcf_parent1 TMP/parent.g.vcf.gz \\
  		 1>&2

if ${keep_gvcf}
then
	mv TMP/child.g.vcf.gz      "${child}.C.g.vcf.gz"
	mv TMP/child.g.vcf.gz.tbi  "${child}.C.g.vcf.gz.tbi"
	mv TMP/parent.g.vcf.gz     "${parent}.P.g.vcf.gz"
	mv TMP/parent.g.vcf.gz.tbi "${parent}.P.g.vcf.gz.tbi"

fi

if ${keep_vcf}
then
	mv TMP/child.vcf.gz        "${child}.C.std.vcf.gz"
	mv TMP/child.vcf.gz.tbi    "${child}.C.std.vcf.gz.tbi"
	mv TMP/parent.vcf.gz       "${parent}.P.std.vcf.gz"
	mv TMP/parent.vcf.gz.tbi   "${parent}.P.std.vcf.gz.tbi"

fi

cat << EOF > versions.yml
"${task.process}":
    deep_variant: "\$(run_deepvariant --version)"
EOF
"""
}
