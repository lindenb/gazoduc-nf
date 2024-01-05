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

include {moduleLoad;getVersionCmd} from '../../modules/utils/functions.nf'


def gazoduc = gazoduc.Gazoduc.getInstance(params)

gazoduc.build("survivor_merge_params","500 1 1 1 0 100").
	menu("Survivor").
        desc("Parameters for survivor (max distance between breakpoints ; Minimum number of supporting caller; Take the type into account;  Take the strands of SVs into accoun; Disabled; Minimum size of SVs to be taken into account).").
        put()


process SURVIVOR_MERGE {
tag "${vcfs.name} contig ${contig}"
afterScript "rm -rf TMP"
memory "20g"
input:
	val(meta)
	path(survivor)
	val(contig) // empty or '*' for any contig
	path(vcfs)
output:
	path("${meta.prefix?:""}survivor.merge.bcf"),emit:vcf
	path("${meta.prefix?:""}survivor.merge.bcf.csi"),emit:index
	path("version.xml"),emit:version
script:
	def survivor_params = meta.survivor_merge_params?:""
	def prefix = meta.prefix?:""
	def viewExtra = contig.equals("*") || contig.isEmpty() ?"":" --regions \"${contig}\" "
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir TMP

i=1
cat "${vcfs}" | while read F
do
	bcftools view ${viewExtra} -O v -o  "TMP/tmp.\${i}.vcf" "\${F}"
	echo "TMP/tmp.\${i}.vcf" >> TMP/sample_files.list
	i=\$((i+1))
done

${survivor.toRealPath()} merge TMP/sample_files.list ${survivor_params} TMP/sample_merged.vcf

	bcftools annotate -x 'INFO/SUPP_VEC' -O u TMP/sample_merged.vcf |\
	bcftools sort --max-mem "${task.memory.giga}G" -T TMP -O u  |\
	bcftools +fill-tags -O u -- -t AN,AC,AF |\
	bcftools view -O b -o "TMP/merged.bcf"

bcftools index TMP/merged.bcf


mv TMP/merged.bcf "${prefix}survivor.merge.bcf"
mv TMP/merged.bcf.csi "${prefix}survivor.merge.bcf.csi"

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">merge VCF file(s) with survivor</entry>
        <entry key="vcfs">${vcfs}</entry>
        <entry key="contig">${contig}</entry>
        <entry key="survivor.version">\$(${survivor}  2>&1 /dev/null | grep '^Version' | cut -d ':' -f 2- )</entry>
        <entry key="survivor.params">${survivor_params} . With:
		* max distance between breakpoints (0-1 percent of length, 1- number of bp) 
		* Minimum number of supporting caller
		* Take the type into account (1==yes, else no)
		* Take the strands of SVs into account (1==yes, else no)
		* Disabled.
		* Minimum size of SVs to be taken into account.</entry>
	<entry key="bcftools.version">${getVersionCmd("bcftools")}</entry>
</properties>
EOF
"""
}
