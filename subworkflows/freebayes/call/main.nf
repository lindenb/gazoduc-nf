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
include {CALL_FREEBAYES as CALL } from '../../../modules/freebayes/call'
include {BCFTOOLS_MERGE         } from '../../../modules/bcftools/merge'
include {BCFTOOLS_CONCAT        } from '../../../modules/bcftools/concat'



workflow FREEBAYES_CALL {
	take:
		meta
		fasta
		fai
		dict
		pedigree
		beds
		bams //[ meta, bam,bai]
	main:
		versions = Channel.empty()
		ch1 = bams
			.map{
				if(it[0].containsKey("batch")) return it;
				return [it[0].plus(batch:it[0].sample),it[1],it[2]];
				}
			.map{[[id:it[0].batch],[it[1],it[2]]]}
			.groupTuple()
			.map{[it[0],it[1].flatten()]}
		
		ch2 = ch1.combine(beds.map{it[1]})

		CALL(
			fasta,
			fai,
			ch2
			)
		versions = versions.mix(CALL.out.versions)


		ch3 = CALL.out.vcf
			.map{[[id:it[3].toRealPath().toString()/* bed */],[it[1],it[2]]]}
			.groupTuple()
			.map{[it[0],it[1].flatten(), [] /* no bed */]}
			.branch{v->
				to_merge:v[1].size()>2 /* bed and it's index */
				other: true
				}


		BCFTOOLS_MERGE(ch3.to_merge)
		versions = versions.mix(BCFTOOLS_MERGE.out.versions)

		without_merge = ch3.other.map{[
			it[0],
			it[1].find{v->v.name.endsWith(".bcf") || v.name.endsWith("*.vcf.gz")}, 
			it[1].find{v->v.name.endsWith(".tbi") || v.name.endsWith("*.csi")}
			]}



		FILTER_AC_GT_0(BCFTOOLS_MERGE.out.vcf.mix(without_merge))
		versions = versions.mix(FILTER_AC_GT_0.out.versions)

	
		BCFTOOLS_CONCAT(
			FILTER_AC_GT_0.out.vcf
				.map{[[id:"freebayes"],[it[1],it[2]]]}
				.groupTuple()
				.map{[it[0],it[1].flatten()]},
			[[id:"nobed"],[]]
			)
		versions = versions.mix(BCFTOOLS_CONCAT.out.versions)


	emit:
		versions
		vcf = BCFTOOLS_CONCAT.out.vcf
	}


process FILTER_AC_GT_0 {
tag "${meta.id}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta),path(vcf),path(vcfidx)
output:
	tuple val(meta),path("*.bcf"),path("*.csi"),emit:vcf
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:vcf.name.md5().substring(0,7)+".ac_gt0"
	// can break java based parser for VCF after BCFTOOLS MERGE: bad GL
	def remove_GL = ((task.ext.no_GL?:false) as boolean)
"""
mkdir -p TMP

bcftools view \\
	-i 'AC[*]>0' \\
	--threads ${task.cpus} \\
	-O u \\
	-o TMP/jeter.bcf \\
	"${vcf}"

if ${remove_GL}
then

	#  FORMAT/GL of freebayes are wrong, crashes htsjdk
	bcftools annotate \\
		--write-index \\
		-x 'FORMAT/GL' \\
		--threads ${task.cpus} \\
		-O b \\
		-o TMP/jeter2.bcf \\
		TMP/jeter.bcf

else

	mv TMP/jeter.bcf TMP/jeter2.bcf

	bcftools index -f \\
		--threads ${task.cpus}  \\
		TMP/jeter2.bcf	

fi




mv TMP/jeter2.bcf ${prefix}.bcf
mv TMP/jeter2.bcf.csi ${prefix}.bcf.csi

cat << EOF > versions.yml
${task.process}:
	bcftools: "\$(bcftools --version | awk '(NR==1){print \$NF}')"
EOF
"""
}
