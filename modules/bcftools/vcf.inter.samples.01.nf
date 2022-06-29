/*

Copyright (c) 2022 Pierre Lindenbaum

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

include {moduleLoad;getKeyValue;assertFileExists;assertNotEmpty} from '../utils/functions.nf'

process VCF_INTER_SAMPLES_01 {
executor "local"
tag "${file(vcf).name} ${file(samples).name}"
afterScript "rm -rf TMP"

input:
	val(meta)
	val(vcf)
	val(samples)
output:
	path("common.samples.txt"),emit:common
	path("vcf_only.txt"),emit:vcf_only
	path("samples_only.txt"),emit:samples_only
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}

mkdir TMP

if [ ! -z "${vcf.endsWith(".list")?"Y":""}" ] ; then

	bcftools concat  --allow-overlaps -O u  --file-list "${vcf}" |\
	bcftools query -l | sort -T TMP | uniq > TMP/jeter.a

else

	bcftools query -l "${vcf}" | sort  -T TMP | uniq > TMP/jeter.a

fi

sort -T TMP "${samples}" | uniq > TMP/jeter.b

comm -12 TMP/jeter.a TMP/jeter.b > common.samples.txt
comm -23 TMP/jeter.a TMP/jeter.b > vcf_only.txt
comm -13 TMP/jeter.a TMP/jeter.b > samples_only.txt

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">get intersection between a VCF and a list of cases/controls</entry>
        <entry key="vcf">${vcf}</entry>
        <entry key="samples">${samples}</entry>
	<entry key="count(vcf.only)">\$( wc -l < vcf_only.txt )</entry>
	<entry key="count(samples.only)">\$( wc -l < samples_only.txt )</entry>
	<entry key="count(common.samples)">\$( wc -l < common.samples.txt )</entry>
	<entry key="bcftools.version">\$( bcftools --version-only )</entry>
</properties>
EOF


# prevent timestamp problems if too fast ?
sleep 5
"""
stub:
"""
touch common.samples.txt vcf_only.txt samples_only.txt
echo "<properties/>" > version.xml
"""
}
