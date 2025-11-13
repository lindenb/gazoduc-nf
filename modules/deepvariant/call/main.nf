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

process DEEP_VARIANT_CALL {
label "process_single"
label "deepvariant"
tag "${meta.id} ${bam.name} ${optional_bed?optional_bed.name:""}"
afterScript "rm -rf TMP  .keras  .parallel"
conda "${moduleDir}/../../../conda/deepvariant.yml"
when:
    task.ext.when == null || task.ext.when
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta ),path(bam),path(bai),path(optional_bed)
output:
	tuple val(meta),path("*.g.vcf.gz"),path("*.g.vcf.gz.tbi"),path(optional_bed),optional:true,emit:gvcf
	tuple val(meta),path("*.std.vcf.gz"),path("*.std.vcf.gz.tbi"),path(optional_bed),optional:true,emit:vcf
	path("versions.yml"),emit:versions
script:

	def sample = meta.id
	def model_type = task.ext.model_type?:""
	if(model_type.isEmpty()) throw new IllegalArgumentException("${task.process} missing ext.model_type (e.g WGS,WES,PACBIO,ONT_R104,HYBRID_PACBIO_ILLUMIN)");
	def prefix = task.ext.prefix?:sample+".deepvariant"
	def keep_gvcf= (task.ext.keep_gvcf?:true).toBoolean()
	def keep_vcf= (task.ext.keep_vcf?:true).toBoolean()
	if(!keep_gvcf && !keep_vcf) throw new IllegalArgumentException("${task.process} discard all");
"""
	hostname 1>&2
	mkdir -p TMP
	mkdir -p TMP/TMP LOGS

	run_deepvariant \\
		--logging_dir LOGS \\
		--model_type ${model_type} \\
		--output_vcf TMP/jeter.vcf.gz \\
		--output_gvcf TMP/jeter.g.vcf.gz \\\
		--num_shards ${task.cpus} \\
		${optional_bed?"--regions ${optional_bed}":""} \\
		--ref "${fasta}" \\
		--reads "${bam}" \\
		--intermediate_results_dir TMP 1>&2

if ${keep_gvcf}
then

		mv TMP/jeter.g.vcf.gz "${prefix}.g.vcf.gz"
		mv TMP/jeter.g.vcf.gz.tbi "${prefix}.g.vcf.gz.tbi"

fi

if ${keep_vcf}
then
	mv TMP/jeter.vcf.gz "${prefix}.std.vcf.gz"
	mv TMP/jeter.vcf.gz.tbi "${prefix}.std.vcf.gz.tbi"
fi

cat << EOF > versions.yml
"${task.process}":
    deep_variant: "\$(run_deepvariant --version)"
EOF
"""
}
