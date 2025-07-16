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

process DEEP_TRIO3 /* 1 child+2 parents */{
label "process_single"
tag "${meta.id} ${C_bam.name} ${optional_bed?optional_bed.name:""}"
afterScript "rm -rf TMP  .keras  .parallel"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta ),path(C_bam),path(C_bai), // CHILD
	 	path(F_bam),path(F_bai),//PARENT 1
		path(M_bam),path(M_bai), //PARENT 2 is optional
		path(optional_bed)
output:
	tuple val(meta),
		path("*.C.g.vcf.gz"),path("*.C.g.vcf.gz.tbi"),
		path("*.F.g.vcf.gz"),path("*.F.g.vcf.gz.tbi"),
		path("*.M.g.vcf.gz"),path("*.M.g.vcf.gz.tbi"),
		path(optional_bed),optional:true,emit:gvcf
	
	tuple val(meta),
		path("*.C.std.vcf.gz"),path("*.C.std.vcf.gz.tbi"),
		path("*.F.std.vcf.gz"),path("*.F.std.vcf.gz.tbi"),
		path("*.M.std.vcf.gz"),path("*.M.std.vcf.gz.tbi"),
		path(optional_bed),optional:true,emit:vcf
		
	path("versions.yml"),emit:versions
script:
	if(!meta.father) throw new IllegalArgumentException("${task.process} missing meta.father")
	if(!meta.mother) throw new IllegalArgumentException("${task.process} missing meta.mother")
	def child = meta.id
	def father = meta.father
	def mother = meta.mother
	def model_type = task.ext.model_type?:(meta.model_type?:"")
	if(model_type.isEmpty()) throw new IllegalArgumentException("${task.process} missing ext.model_type (e.g WGS,WES,PACBIO,ONT_R104,HYBRID_PACBIO_ILLUMINA)");
	def keep_gvcf= (task.ext.keep_gvcf?:true) as boolean
	def keep_vcf= (task.ext.keep_vcf?:true) as boolean
"""
	hostname 1>&2
	mkdir -p TMP/TMP TMP/LOGS
	
	run_deeptrio \\
		--logging_dir TMP/LOGS \\
		--model_type ${model_type} \\
		--ref "${fasta}" \\
		--reads_child "${C_bam}" \\
		--reads_parent1 "${F_bam}" \\
		--reads_parent2 "${M_bam}" \\
		--sample_name_child "${child}" \\
	    --sample_name_parent1 "${father}" \\
	    --sample_name_parent2 "${mother}" \\
		--num_shards ${task.cpus} \\
		--intermediate_results_dir TMP/TMP \\
		${optional_bed?"--regions ${optional_bed}":""} \\
		--output_vcf_child TMP/child.vcf.gz \\
  		--output_vcf_parent1 TMP/father.vcf.gz \\
  		--output_vcf_parent2 TMP/mother.vcf.gz \\
  		--output_gvcf_child  TMP/child.g.vcf.gz \\
  		--output_gvcf_parent1 TMP/father.g.vcf.gz \\
  		--output_gvcf_parent2 TMP/mother.g.vcf.gz \\
		1>&2


	find TMP -type f 1>&2

if ${keep_gvcf}
then
	mv TMP/child.g.vcf.gz      "${child}.C.g.vcf.gz"
	mv TMP/child.g.vcf.gz.tbi  "${child}.C.g.vcf.gz.tbi"
	mv TMP/father.g.vcf.gz     "${father}.F.g.vcf.gz"
	mv TMP/father.g.vcf.gz.tbi "${father}.F.g.vcf.gz.tbi"
	mv TMP/mother.g.vcf.gz     "${mother}.M.g.vcf.gz"
	mv TMP/mother.g.vcf.gz.tbi "${mother}.M.g.vcf.gz.tbi"
fi

if ${keep_vcf}
then
	mv TMP/child.vcf.gz        "${child}.C.std.vcf.gz"
	mv TMP/child.vcf.gz.tbi    "${child}.C.std.vcf.gz.tbi"
	mv TMP/father.vcf.gz       "${father}.F.std.vcf.gz"
	mv TMP/father.vcf.gz.tbi   "${father}.F.std.vcf.gz.tbi"
	mv TMP/mother.vcf.gz       "${mother}.M.std.vcf.gz"
	mv TMP/mother.vcf.gz.tbi   "${mother}.M.std.vcf.gz.tbi"
fi

cat << EOF > versions.yml
"${task.process}":
    deep_trio: "\$(run_deeptrio --version 2> /dev/null  | awk '(NR==1) {print \$NF;}')"
EOF
"""
}
