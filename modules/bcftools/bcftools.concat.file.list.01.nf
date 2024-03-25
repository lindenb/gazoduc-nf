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
include {moduleLoad} from '../../modules/utils/functions.nf'


process BCFTOOL_CONCAT_FILE_LIST_01 {
tag "${vcfs.name}"
afterScript "rm -rf TMP"
cpus 1
input:
	val(meta)
        path(vcfs)
	path(bed) // or NO_FILE
output:
        path("${params.prefix?:""}concatenated.bcf"),emit:vcf
        path("${params.prefix?:""}concatenated.bcf.csi"),emit:index
	path("version.xml"),emit:version
script:
	def rgn = bed.name.equals("NO_FILE")?"":"--regions-file \"${bed}\""
	def concat_args = " ${rgn} --no-version --allow-overlaps --remove-duplicates -O b "
	def limit = 200
"""
	hostname 1>&2
	${moduleLoad("bcftools")}
	mkdir -p TMP
	set -x
	if test \$(wc -l < "${vcfs}") -gt ${limit}
	then
		cat "${vcfs}" | paste -d, - - - - - - - - - - - | while read F
		do
			echo "# \${F}" 1>&2
			if ! test -f TMP/jeter1.bcf
			then
				bcftools concat --threads ${task.cpus} ${concat_args} -o TMP/jeter2.bcf  \${F//,/ }				
			else
				bcftools concat --threads ${task.cpus} ${concat_args} -o TMP/jeter2.bcf  TMP/jeter1.bcf \${F//,/ }
			fi
			mv -v TMP/jeter2.bcf TMP/jeter1.bcf
			bcftools index --force --threads ${task.cpus} TMP/jeter1.bcf
		done

		mv -v  TMP/jeter1.bcf "${params.prefix?:""}concatenated.bcf"
		mv -v TMP/jeter1.bcf.csi "${params.prefix?:""}concatenated.bcf.csi"
	else

		bcftools concat --threads ${task.cpus} ${concat_args} -o "TMP/${params.prefix?:""}concatenated.bcf" --file-list "${vcfs}"

		bcftools index --threads ${task.cpus}  "TMP/${params.prefix?:""}concatenated.bcf"

		mv -v TMP/${params.prefix?:""}concatenated.* ./
	fi


	##################################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">concat vcf(s) using bcftools</entry>
		<entry key="vcfs">${vcfs}</entry>
		<entry key="count(vcfs)">\$(wc -l < ${vcfs})</entry>
		<entry key="bcftools">\$( bcftools --version-only)</entry>
	</properties>
	EOF
"""
}
