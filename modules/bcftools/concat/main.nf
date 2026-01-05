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


process BCFTOOLS_CONCAT {
label "process_short"
tag "${meta.id?:""}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.02.yml"
input:
	tuple val(meta ),path("VCFS/*"),path(optional_bed)
output:
    tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),path(optional_bed),emit:vcf
    tuple val(meta),path("*.md5"),optional:true,emit:md5
	path("versions.yml"),emit:versions
script:
	def args1 = task.ext.args1?:""
	def args2 = task.ext.args2?:"--no-version --allow-overlaps --remove-duplicates "
	def args3 = optional_bed?"--regions-file \"${optional_bed}\"":""
	def limit = task.ext.limit?:10
	def prefix = task.ext.prefix?:(meta.id?:"variants")+".\${MD5}.concat"
	
	def by_contig = (task.ext.by_chromosome?:false).toBoolean()
	def with_md5 = (task.ext.with_md5?:true).toBoolean()
	log.warn("BCFTOOLS_CONCAT deprecated")
"""	
	hostname 1>&2
	mkdir -p TMP
	find VCFS/ \\( -name "*.vcf.gz" -o -name "*.bcf" \\) | sort -V -T TMP > TMP/jeter.list
	MD5=\$(cat TMP/jeter.list ${optional_bed?optional_bed:""} | md5sum | cut -d ' ' -f1)

	set -x

	if ${by_contig}
	then
		cat TMP/jeter.list | while read F
		do
			bcftools index --threads ${task.cpus}  -s "\${F}" |\\
				awk -F '\t' -vF="\$F" '(printf("%s\t0\t%s\t%s\\n",\$1,\$2,F);}' |\\
				sort -T TMP -k1,1 -k2,2n >> TMP/contigs.bed
		done

		# each distinct contig
		cut -f1 TMP/contigs.bed | sort | uniq | while read CONTIG
		do
			# build a BED for this CONTIG
			awk -F '\t' -vC="\${CONTIG}" '(\$1==C)' TMP/contigs.bed |\\
				cut -f1,2,3 |\\
				sort -T TMP -k1,1 -k2,2n |\\
				bedtools merge > "TMP/\${CONTIG}.ctg.bed"

			# intersect with optional bed, if any
			if ${optional_bed?true:false}
			then
				bedtools intersect -a "${optional_bed}" -b "TMP/\${CONTIG}.ctg.bed" |\\
				sort -T TMP -k1,1 -k2,2n |\\
				bedtools merge > "TMP/\${CONTIG}.ctg.bed.new"

				mv "TMP/\${CONTIG}.ctg.bed.new"  "TMP/\${CONTIG}.ctg.bed"
			fi

			# no overlap with optional bed
			if test ! -s "TMP/\${CONTIG}.ctg.bed"
			then
				continue
			fi

			#vcf overlaping this contigs
			awk -F '\t' -vC="\${CONTIG}" '(\$1==C)' TMP/contigs.bed |\\
				cut -f4 | sort | uniq > TMP/jeter2.vcf.list

			test -s TMP/jeter2.vcf.list

			# too many files or not ?
			if test  `wc -l < TMP/jeter2.vcf.list` -le ${limit}
					then

						bcftools concat \\
							--threads ${task.cpus} \\
							${args1} \\
							${args2} \\
							--regions-file "TMP/\${CONTIG}.ctg.bed" \\
							-O z9 \\
							--file-list TMP/jeter2.vcf.list \\
							-o "TMP/jeter.vcf.gz" 

					else
						
						SQRT=`awk 'END{X=NR;if(X<10){print(X);} else {z=sqrt(X); print (z==int(z)?z:int(z)+1);}}' TMP/jeter2.vcf.list`
						rm -fv TMP/chunck*
						rm -fv TMP/jeter3.list
						split -a 9 --additional-suffix=.list --lines=\${SQRT} TMP/jeter2.vcf.list TMP/chunck.


						find TMP/ -type f -name "chunck*.list" | while read F
						do
							bcftools concat \\
								--write-index \\
								--threads ${task.cpus} \\
								${args1} \\
								${args2} \\
								--regions-file "TMP/\${CONTIG}.ctg.bed" \\
								-O b \\
								--file-list "\${F}" \\
								-o "\${F}.delete.me.bcf" 
							echo "\${F}.delete.me.bcf" >> TMP/jeter3.list
						done

						bcftools concat \\
							--threads ${task.cpus} \\
							${args1} \\
							${args2} \\
							--regions-file "TMP/\${CONTIG}.ctg.bed" \\
							-O z9 \\
							--file-list TMP/jeter3.list \\
							-o "TMP/jeter.vcf.gz" 

						rm -vf "TMP/*.delete.me.bcf"
						rm -vf "TMP/*.delete.me.bcf.csi"
					fi

					# default write index is CSI, not TBI for vcf.gz
				
					bcftools index -f -t --threads ${task.cpus} TMP/jeter.vcf.gz
					

					mv -v TMP/jeter.vcf.gz "${prefix}.\${CONTIG}.vcf.gz"
					mv -v TMP/jeter.vcf.gz.tbi "${prefix}.\${CONTIG}.vcf.gz.tbi"

					# Generate MD5 if needed
					if ${with_md5}
					then 
						md5sum "${prefix}.\${CONTIG}.vcf.gz" > "${prefix}.\${CONTIG}.vcf.gz.md5"
					fi



		done
		# end of loop over each chromosome


	else
		# do NOT group by chromosome
		
		if test  `wc -l < TMP/jeter.list` -eq 1
		then

			bcftools view \\
				--threads ${task.cpus} \\
				${args3} \\
				-O z9 \\
				-o "TMP/jeter.vcf.gz" \\
				`cat TMP/jeter.list`

		elif test  `wc -l < TMP/jeter.list` -le ${limit}
		then

			bcftools concat \\
				--threads ${task.cpus} \\
				${args1} \\
				${args2} \\
				${args3} \\
				-O z9 \\
				--file-list TMP/jeter.list \\
				-o "TMP/jeter.vcf.gz"

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
				--threads ${task.cpus} \\
				${args1} \\
				${args2} \\
				${args3} \\
				-O z9 \\
				--file-list TMP/jeter2.list \\
				-o "TMP/jeter.vcf.gz"

		fi

		# default write index is CSI, not TBI for vcf.gz
		bcftools index -f -t --threads ${task.cpus} TMP/jeter.vcf.gz

		mv TMP/jeter.vcf.gz ${prefix}.vcf.gz
		mv TMP/jeter.vcf.gz.tbi ${prefix}.vcf.gz.tbi


		# Generate MD5 if needed
		if ${with_md5}
		then 
			md5sum ${prefix}.vcf.gz > ${prefix}.vcf.gz.md5
		fi

	fi



cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""

stub:

"""
find VCFS/ \\( -name "*.vcf.gz" -o -name "*.bcf" \\) | sort -V -T . 1>&2
touch versions.yml ${meta.id}.vcf.gz  ${meta.id}.vcf.gz.tbi   ${meta.id}.vcf.gz.md5
"""
}
