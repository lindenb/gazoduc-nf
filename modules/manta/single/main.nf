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

process MANTA_SINGLE {
    tag "${meta.id?:""}"
    label "process_single"
    label "manta"
    afterScript "rm -rf TMP"
    conda "${moduleDir}/../../../conda/manta.yml"
    input:
		tuple val(meta1),path(fasta)
		tuple val(meta2),path(fai)
		tuple val(meta3),path(dict)
        tuple val(meta4),path(optional_bed)
        tuple val(meta),path(bam),path(bai)
    output:
        tuple val(meta),path("*.diploidSV.vcf.gz"),path("*.diploidSV.vcf.gz.tbi"),emit:diploidSV
        tuple val(meta),path("*.candidateSmallIndels.vcf.gz"),path("*.candidateSmallIndels.vcf.gz.tbi"),emit:candidateSmallIndels
        tuple val(meta),path("*.candidateSV.vcf.gz"),path("*.candidateSV.vcf.gz.tbi"),emit:candidateSV
    	path("versions.yml"),emit:versions
    script:
        def has_bed = optional_bed?true:false
        def args1 = task.ext.args1?:""
        def args2 = task.ext.args2?:"--quiet"
        def prefix = task.ext.prefix?:meta.id
	"""
	hostname 1>&2
        mkdir -p TMP
        set -x
 	echo "\${PATH}" 1>&2

	if ${has_bed}
    then
		${has_bed && optional_bed.name.endsWith(".gz")?"gunzip -c":"cat"} ${optional_bed} |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n > manta.bed
	    bgzip manta.bed
	    tabix -fp bed manta.bed.gz
	fi


	configManta.py \\
		--bam ${bam} \\
        ${args1} \\
		--referenceFasta "${fasta}" \\
		${has_bed?"--callRegions manta.bed.gz":""} \\
		--runDir "TMP"


	./TMP/runWorkflow.py ${args2} -m local -j ${task.cpus}
	rm -rf TMP/workspace
	
	for X in diploidSV candidateSmallIndels candidateSV
	do
		bcftools view --threads ${task.cpus} -O z -o "${prefix}.\${X}.vcf.gz" "./TMP/results/variants/\${X}.vcf.gz"
		bcftools index --threads ${task.cpus}  -f -t "${prefix}.\${X}.vcf.gz"
	done


cat << EOF > versions.yml
"${task.process}":
	manta: "\$(configManta.py --version)"
EOF
	"""


stub:
	def prefix = "${meta.id}"
"""
touch versions.yml

for X in diploidSV candidateSmallIndels candidateSV
do
	touch "${prefix}.\${X}.vcf.gz"
	touch "${prefix}.\${X}.vcf.gz.tbi"
done

"""
}
