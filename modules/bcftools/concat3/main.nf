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
include {parseBoolean } from '../../utils/functions.nf'

/**
 *
 *  simple bcftools concat, producing only one vcf, no bed
 *
 */
process BCFTOOLS_CONCAT {
label "process_short"
tag "${meta.id?:""}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.02.yml"
input:
	tuple val(meta ),path("VCFS/*")
output:
	tuple val(meta),path("*.vcf.gz"),path("*.tbi"),emit:vcf
	tuple val(meta),path("*.md5"),optional:true,emit:md5
	path("versions.yml"),emit:versions
script:
	def args1 = task.ext.args1?:""
	def limit = (task.ext.limit?:100) as int
	def args2 = task.ext.args2?:"--no-version --allow-overlaps --remove-duplicates "
	def prefix = task.ext.prefix?:"${meta.id}.concat"	
	def with_md5 =  parseBoolean(task.ext.with_md5?:true)
"""	
	hostname 1>&2
	mkdir -p TMP
	
	find VCFS/ \\( -name "*.vcf.gz" -o -name "*.bcf" \\) | sort -V -T TMP > TMP/jeter.list

	my_concat () {
		bcftools concat \\
			--threads ${task.cpus} \\
			--write-index \\
			${args1} \\
			${args2} \\
			--file-list "\${1}" \\
			-O "\${2}" \\
			-o "\${3}"
		}


	set -x
	if test  `wc -l < TMP/jeter.list` -le ${limit}
	then
		my_concat TMP/jeter.list z9 TMP/jeter.vcf.gz
	else

		SQRT=`awk 'END{X=NR;if(X<10){print(X);} else {z=sqrt(X); print (z==int(z)?z:int(z)+1);}}' TMP/jeter.list`
		split -a 9 --additional-suffix=.list --lines=\${SQRT} TMP/jeter.list TMP/chunck.


		find TMP/ -type f -name "chunck*.list" | while read F
		do
			my_concat "\${F}" b "\${F}.bcf"
			
			echo "\${F}.bcf" >> TMP/jeter2.list
		done

		test -s TMP/jeter2.list
		
		my_concat TMP/jeter2.list z9 TMP/jeter.vcf.gz

	fi


	# default write index is CSI, not TBI for vcf.gz
	bcftools index -f -t --threads ${task.cpus} TMP/jeter.vcf.gz

	mv TMP/jeter.vcf.gz ${prefix}.vcf.gz
	mv TMP/jeter.vcf.gz.tbi ${prefix}.vcf.gz.tbi


	# Generate MD5 if needed
	if ${with_md5}
	then 
		md5sum ${prefix}.vcf.gz > ${prefix}.vcf.gz.md5
	fi


cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""

stub:

"""
find VCFS/ \\( -name "*.vcf.gz" -o -name "*.bcf" \\) | sort -V -T . 1>&2
touch versions.yml ${meta.id}.vcf.gz  ${meta.id}.vcf.gz.tbi   ${meta.id}.vcf.gz.md5
"""
}
