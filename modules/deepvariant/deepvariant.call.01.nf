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

include {moduleLoad} from '../utils/functions.nf'

process DEEP_VARIANT_CALL_01 {
tag "${row.sample}/${file(row.bed).name}/${file(row.bam).name}"
afterScript "rm -rf TMP  .keras  .parallel"
input:
	val(meta)
	val(row)
output:
	tuple val(row),path("${row.sample}.bcf"),path("${row.sample}.g.vcf.gz"),emit:output
	path("version.xml"),emit:version
script:
	def img = params.deepvariant.singularity_image
	def sample = row.sample
	def new_sample = row.new_sample?:sample
	def bam = file(row.bam)
	def bed = file(row.bed)
	def genome = params.genomes[row.genomeId]
	def ref = file(genome.fasta)
	def model_type = row.model_type?:params.deepvariant.model_type
"""
	hostname 1>&2
	${moduleLoad("bcftools")}
	mkdir TMP
	
	singularity exec\
		--home \${PWD} \
		--bind ${bam.getParent()}:/data1 \
		--bind ${bed.getParent()}:/beddir \
		--bind `dirname ${ref.toRealPath()}`:/ref \
		--bind \${PWD}/TMP:/outdir \
		${img} \
		/opt/deepvariant/bin/run_deepvariant \
		--model_type ${model_type} \
		--output_vcf /outdir/jeter.vcf.gz \
		--output_gvcf /outdir/jeter.g.vcf.gz \
		--num_shards ${task.cpus} \
		--regions /beddir/${bed.name} \
		--ref /ref/${ref.name} \
		--reads /data1/${bam.name} \
		--intermediate_results_dir  /outdir/TMP 1>&2



	mv TMP/jeter.g.vcf.gz "${sample}.g.vcf.gz"
	mv TMP/jeter.g.vcf.gz.tbi "${sample}.g.vcf.gz.tbi"

	bcftools view -O b -o "${sample}.bcf" TMP/jeter.vcf.gz
	bcftools index "${sample}.bcf"

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">call bam with deepvariant</entry>
        <entry key="sample">${sample}</entry>
        <entry key="bam">${bam}</entry>
        <entry key="model_type">${model_type}</entry>
        <entry key="deepvariant.version">\$(singularity exec ${img} /opt/deepvariant/bin/run_deepvariant --version) </entry>
</properties>
EOF
"""
}
