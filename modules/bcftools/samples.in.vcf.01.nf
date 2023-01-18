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

include {getVersionCmd;moduleLoad;getKeyValue;assertFileExists;assertNotEmpty} from '../utils/functions.nf'

/**
extract samples from VCF using bcftools query
**/
process SAMPLES_IN_VCF_01 {
executor "local"
tag "${file(vcf).name}"

input:
	val(meta)
	val(vcf)
output:
	path("samples.txt"),emit:samples
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}
set -x
#set -o pipefail

if ${vcf.endsWith(".list")} ; then

	bcftools concat  --allow-overlaps -O u  --file-list "${vcf}" |\
	bcftools query -l | sort | uniq > jeter.a

else

	bcftools query -l "${vcf}" | sort | uniq > jeter.a

fi

mv -v jeter.a "samples.txt"

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">extract samples from vcf</entry>
        <entry key="vcf">${vcf}</entry>
        <entry key="samples.count">\$(wc -l < samples.txt )</entry>
	<entry key="versions">${getVersionCmd("awk bcftools")}</entry>
</properties>
EOF
"""
}

