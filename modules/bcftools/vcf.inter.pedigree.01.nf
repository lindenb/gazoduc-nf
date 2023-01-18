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

include {moduleLoad;getKeyValue;assertFileExists;assertNotEmpty} from '../utils/functions.nf'

process VCF_INTER_PED_01 {
tag "${file(vcf).name} ${file(pedigree).name}"

afterScript "rm -f jeter.a jeter.b jeter.samples"

input:
	val(meta)
	val(vcf)
	val(pedigree)
output:
	path("${prefix}samples.common.txt"),emit:samples
	path("${prefix}pedigree.common.ped"),emit:pedigree
	path("${prefix}invcf.notinped.txt"),emit:notinped
	path("${prefix}inped.notinvcf.txt"),emit:notnvcf
	path("version.xml"),emit:version
script:
	prefix = getKeyValue(meta,"prefix","")
	ped_type = assertNotEmpty(getKeyValue(meta,"pedigree_type",""),"pedigree_type must be defined")
	assertFileExist(vcf,"vcf must be defined")
	assertFileExist(vcf,"pedigree must be defined")
"""
hostname 1>&2
${moduleLoad("bcftools")}
set -x

if [ ! -z "${vcf.endsWith(".list")?"Y":""}" ] ; then

	bcftools concat  --allow-overlaps -O u  --file-list "${vcf}" |\
	bcftools query -l | sort | uniq > jeter.a

else

	bcftools query -l "${vcf}" | sort | uniq > jeter.a

fi

grep -v "^fid\tiid" "${pedigree}" | cut -f2 | sort | uniq > jeter.b
comm -12 jeter.a jeter.b > jeter.samples
# test will fail if there is no common samples between vcf and the pedigree
test -s jeter.samples


comm -23 jeter.a jeter.b > "${prefix}invcf.notinped.txt"
comm -13 jeter.a jeter.b > "${prefix}inped.notinvcf.txt"


# save header
if [ "${ped_type}" == "rvtests" ] ;
	grep -m1 "^fid\tiid" "${pedigree}" > jeter2.ped
else
	touch jeter2.ped
fi

# join with pedigree/samples, discard samples without phenotype
grep -v "^fid\tiid" "${pedigree}" |\
	awk '(\$6=="1" || \$6=="2")' |\
	sort -t '\t' -T . -k2,2 |\
	join -t '\t' -1 1 -2 2 -o '2.1,2.2,2.3,2.4,2.5,2.6' jeter.samples - >> jeter2.ped





mv -v jeter.samples "${prefix}samples.common.txt"
mv -v jeter2.ped "${prefix}pedigree.common.ped"

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">get intersection between a VCF and a pedigree</entry>
        <entry key="pedigree.type">${ped_type}</entry>
        <entry key="vcf">${vcf}</entry>
        <entry key="pedigree">${pedigree}</entry>
        <entry key="samples.not.in.vcf">\$(paste -d ' ' -s < "${prefix}inped.notinvcf.txt" )</entry>
        <entry key="samples.not.in.ped">\$(paste -d ' ' -s < "${prefix}invcf.notinped.txt" )</entry>
	<entry key="bcftools.version">\$( bcftools --version-only)</entry>
</properties>
EOF

"""
stub:
"""
touch "${prefix}samples.common.txt"
touch "${prefix}pedigree.common.ped"
touch "${prefix}invcf.notinped.txt"
touch "${prefix}inped.notinvcf.txt"
echo "<properties/>" > version.xml
"""
}

