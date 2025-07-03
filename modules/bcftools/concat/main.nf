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


process BCFTOOL_CONCAT {
label "process_single"
tag "${meta.id?:""}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	tuple val(meta ),path("VCFS/*") 
	tuple val(meta1),path(optional_bed)
output:
    tuple val(meta),path("*.{bcf,vcf.gz}"),path("*.{csi,tbi}"),emit:vcf
	path("versions.yml"),emit:versions
script:
	def args1 = task.ext.args1?:""
	def args2 = task.ext.args2?:"--no-version --allow-overlaps --remove-duplicates "
	def args3 = optional_bed?"--regions-file \"${optional_bed}\"":""
	def limit = task.ext.limit?:10
	def prefix = task.ext.prefix?:meta.id+".concat"
"""	
	hostname 1>&2
	mkdir -p TMP
	find VCFS/ -name "*.vcf.gz" -o -name "*.bcf" > TMP/jeter.list


	if test  `wc -l < TMP/jeter.list` -le ${limit}
	then

		bcftools concat \\
			--write-index ${args} \\
			--threads ${task.cpus} \\
			${args1} \\
			${args2} \\
			${args3} \\
			--O b9 \\
			--file-list TMP/jeter.list \\
			-o "TMP/jeter.bcf" 

	else
	
		SQRT=`awk 'END{X=NR;if(X<10){print(X);} else {z=sqrt(X); print (z==int(z)?z:int(z)+1);}}' TMP/jeter.list`
		split -a 9 --additional-suffix=.list --lines=\${SQRT} TMP/jeter.list TMP/chunck.


		find TMP/ -type f -name "chunck*.list" | while read F
		do
			bcftools concat \\
				--write-index \\
				--threads ${task.cpus} \\
				${args1} \\
				${args2} \\
				${args3} \\
				-O b \\
				--file-list "\${F}" \\
				-o "\${F}.bcf" 
			echo "\${F}.bcf" >> TMP/jeter2.list
		done

		bcftools concat \\
			--write-index \\
			--threads ${task.cpus} \\
			${args1} \\
			${args2} \\
			${args3} \\
			--O b9 \\
			--file-list TMP/jeter2.list \\
			-o "TMP/jeter.bcf" 

	fi

	mv TMP/jeter.bcf ${prefix}.bcf
	mv TMP/jeter.bcf.csi ${prefix}.bcf.csi

cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
}
