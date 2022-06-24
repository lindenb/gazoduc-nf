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
include {moduleLoad;getKeyValue} from '../../modules/utils/functions.nf'

process BCFTOOLS_MERGE_BED_01 {
tag "${file(bed).name} N=${L.size()}"
afterScript "rm -f jeter3.list jeter2.list jeter.list"
cpus 1
input:
	val(meta)
	tuple val(bed),val(L)
output:
       	path("genotyped.bcf"),emit:bed
        path("genotyped.bcf.csi"),emit:index
	path("version.xml"),emit:version
script:
	def extraMerge = getKeyValue(meta,"extraMerge","")
	def extraView = getKeyValue(meta,"extraView","")

if(L.size()==1) 
"""
	hostname 1>&2
	${moduleLoad("bcftools")}

	bcftools view --threads ${task.cpus} ${extraView} --regions-file "${bed}" -O b -o genotyped.bcf "${L[0]}"
	bcftools index --threads ${task.cpus} "genotyped.bcf"


	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">merge invoked, but there was only one vcf: just copy.</entry>
		<entry key="bed">${bed}</entry>
		<entry key="vcf-count">${L.size()}</entry>
		<entry key="extra view params">${extraView}</entry>
		<entry key="bcftools">\$( bcftools --version-only)</entry>
	</properties>
	EOF
"""
else
"""
hostname 1>&2
${moduleLoad("bcftools")}
set -o pipefail

cat << EOF > jeter.list
${L.join("\n")}
EOF

#
# use the first sample of each VCF to be sure that samples will be ordered the same way 
# for any downstream bcftools concat
#
cat jeter.list | while read V
do
	bcftools query -l "\$V" | awk -vF=\$V 'BEGIN{FOUND=0;}(NR==1) {printf("%s,%s\\n",\$1,F);FOUND=1;} END{if(FOUND==0) {printf("NO_GT,%s\\n",F);} }' >> jeter2.list
done

LC_ALL=C sort -t, -T . -k1,1 jeter2.list | cut -d, -f2 > jeter3.list

bcftools merge --threads ${task.cpus} ${extraMerge} ${extraView} --regions-file "${bed}"  --no-version -O b -o "genotyped.bcf" --file-list jeter3.list
bcftools index --threads ${task.cpus} "genotyped.bcf"


	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">merge VCFs</entry>
		<entry key="bed">${bed}</entry>
		<entry key="vcf-count">${L.size()}</entry>
		<entry key="extra view params">${extraView}</entry>
		<entry key="extra merge params">${extraMerge}</entry>
		<entry key="bcftools">\$( bcftools --version-only)</entry>
	</properties>
	EOF



"""
stub:
"""
touch genotyped.bcf genotyped.bcf
echo "<properties/>" > version.xml
"""
}
