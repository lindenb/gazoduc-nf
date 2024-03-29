/*

Copyright (c) 2023 Pierre Lindenbaum

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
include {moduleLoad;removeCommonSuffixes} from '../../modules/utils/functions.nf'

process BCFTOOL_STATS_01 {
tag "${vcf.name}"
afterScript "rm -rf TMP"
cpus 1
input:
        val(meta)
	val(reference)
        path(vcf)
output:
        tuple path(vcf),path("${meta.prefix?:""}${removeCommonSuffixes(vcf.name)}.zip"),emit:output
	path("version.xml"),emit:version
script:
	def name = (meta.prefix?:"") + removeCommonSuffixes(vcf.name)
	def args =  "" + 
		(meta.samples_file?"--samples-file '"+meta.samples_file+"'":"") + 
		(meta.samples?"--samples '"+ meta.samples+"'":"")
"""
	hostname 1>&2
	${moduleLoad("bcftools")}
	set -o pipefail
	mkdir TMP

	bcftools stats --threads ${task.cpus} -F ${reference} ${args} ${vcf.toRealPath()}  > 'TMP/${name}.stats.txt'

	# ugly !!!!!!!!!!!!!!
	cat <<- EOF | bash -l -e
	conda activate matplotlib
	plot-vcfstats -p TMP TMP/${name}.stats.txt
	conda deactivate
	EOF

	gzip --best 'TMP/${name}.stats.txt'

	rm -f TMP/*.aux
	mv TMP "${name}"
	
	zip -r "${name}.zip" "${name}/"
	

	rm -rf "${name}"

	##################################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">bcftools stats</entry>
		<entry key="vcf">${vcf.name}</entry>
		<entry key="args"><code><![CDATA[${args}]]></code></entry>
		<entry key="bcftools">\$( bcftools --version-only)</entry>
	</properties>
	EOF
"""
}
