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
process BCFTOOLS_MERGE {
label "process_short"
tag "${meta.id?:""}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta ),path("VCFS/*"),path(optional_bed)
output:
        tuple val(meta),path("*.bcf"),path("*.bcf.csi"),emit:vcf
	path("versions.yml"),emit:versions
script:
        def args1  = task.ext.args1?:""
        def prefix = task.ext.prefix?:"\${MD5}"+(meta.id?"."+meta.id.md5().substring(0,7):"") + (optional_bed?"."+optional_bed.baseName:"")
"""
mkdir -p TMP
find VCFS/ \\( -name "*.vcf.gz" -o -name "*.bcf" \\)  > TMP/jeter.list

#
# use the first sample of each VCF to be sure that samples will be ordered the same way 
# for any downstream bcftools concat
#
cat TMP/jeter.list | while read V
do
	bcftools query -l "\$V" |\\
		awk -vF=\$V 'BEGIN{FOUND=0;}(NR==1) {printf("%s,%s\\n",\$1,F);FOUND=1;} END{if(FOUND==0) {printf("NO_GT,%s\\n",F);} }' >> TMP/jeter2.list
done

LC_ALL=C sort -t, -T TMP -k1,1 TMP/jeter2.list | cut -d, -f2 > TMP/jeter3.list
MD5=`cat TMP/jeter3.list | md5sum | cut -d ' ' -f1`

if [[ \$(wc -l < TMP/jeter.list) -eq 1 ]]
then

bcftools view \\
	--threads ${task.cpus} \\
	${optional_bed?"--regions-file \"${optional_bed}\"":""} \\
	-O u \\
	-o TMP/jeter2.bcf \\
	`cat TMP/jeter.list`


else

bcftools merge \\
	--threads ${task.cpus} \\
	${args1} \\
	${optional_bed?"--regions-file \"${optional_bed}\"":""} \\
	--no-version \\
	-O u \\
	-o "TMP/jeter2.bcf" \\
	--file-list TMP/jeter3.list
fi

bcftools  +fill-tags \\
	--threads ${task.cpus} \\
	-O b  \\
	-o TMP/jeter.bcf \\
	TMP/jeter2.bcf -- -t AN,AC,AF,AC_Hom,AC_Het,AC_Hemi,NS

bcftools index  -f --threads ${task.cpus}  TMP/jeter.bcf

mv TMP/jeter.bcf     ${prefix}.bcf
mv TMP/jeter.bcf.csi ${prefix}.bcf.csi

cat << END_VERSIONS > versions.yml
${task.process}:
    bcftools: \$(bcftools version | awk '(NR==1)  {print \$NF}')
END_VERSIONS
"""
}
