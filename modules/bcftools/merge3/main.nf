/*

Copyright (c) 2026 Pierre Lindenbaum

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


/*
 *
 * ultra simple bcftools merge, no BED
 *
 */
process BCFTOOLS_MERGE {
label "process_single"
tag "${meta.id?:""}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta ),path("VCFS/*")
output:
        tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
	path("versions.yml"),emit:versions
script:
        def args1  = task.ext.args1?:""
        def args2  = task.ext.args2?:""
        def args3  = task.ext.args3?:""
        def prefix = task.ext.prefix?:"${meta.id}"
		def cmd = task.ext.cmd?:"view"
		def tags = task.ext.tags?:"AN,AC,AF,AC_Hom,AC_Het,AC_Hemi,NS"
		def limit = (task.ext.limit?:100) as int
"""
mkdir -p TMP
find VCFS/ \\( -name "*.vcf.gz" -o -name "*.bcf" \\)  > TMP/jeter.list


my_merge () {
		bcftools merge \\
			--threads ${task.cpus} \\
			--write-index \\
			${args1} \\
			--file-list "\${1}" \\
			-O "\${2}" \\
			-o "\${3}"
		}

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

rm -f TMP/jeter2.list

if [[ \$(wc -l < TMP/jeter3.list) -eq 1 ]]
then

	bcftools view \\
		--threads ${task.cpus} \\
		-O u \\
		-o TMP/jeter2.bcf \\
		`cat TMP/jeter.list`

elif test  `wc -l < TMP/jeter3.list` -le ${limit}
then

	my_merge TMP/jeter3.list u TMP/jeter2.bcf

else
		
		SQRT=`awk 'END{X=NR;if(X<10){print(X);} else {z=sqrt(X); print (z==int(z)?z:int(z)+1);}}' TMP/jeter3.list`
		split -a 9 --additional-suffix=.list --lines=\${SQRT} TMP/jeter3.list TMP/chunck.


		find TMP/ -type f -name "chunck*.list" | while read F
		do
			my_merge "\${F}" b "\${F}.bcf"
			
			echo "\${F}.bcf" >> TMP/jeter2.list
		done

		test -s TMP/jeter2.list
		
		my_merge TMP/jeter2.list u TMP/jeter2.bcf

fi

# give a chance to filter-out things, or bcftools annotate
bcftools ${cmd} ${args2} -O u TMP/jeter2.bcf |\\
bcftools  +fill-tags \\
	${args3} \\
	--threads ${task.cpus} \\
	-O z9  \\
	-o TMP/jeter.vcf.gz \\
	-- -t ${tags}

bcftools index  -f -t --threads ${task.cpus}  TMP/jeter.vcf.gz

mv TMP/jeter.vcf.gz  ${prefix}.vcf.gz
mv TMP/jeter.vcf.gz.tbi  ${prefix}.vcf.gz.tbi

cat << END_VERSIONS > versions.yml
${task.process}:
    bcftools: \$(bcftools version | awk '(NR==1)  {print \$NF}')
END_VERSIONS
"""

stub:
        def prefix = task.ext.prefix?:"${meta.id}"
 // "f"+(meta.id?"."+meta.id.md5():"") + (optional_bed?"."+optional_bed.baseName:"")
"""
find VCFS/ \\( -name "*.vcf.gz" -o -name "*.bcf"  \\) 
touch versions.yml ${prefix}.vcf.gz ${prefix}.vcf.gz.tbi
"""
}
