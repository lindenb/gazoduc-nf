/*

Copyright (c) 2025 Pierre Lindenbaum

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
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta ),path("VCFS/*") 
	tuple val(meta1),path(optional_bed)
output:
    tuple val(meta),path("*.{bcf,vcf.gz}"),path("*.{csi,tbi}"),optional:true,emit:vcf
	path("versions.yml"),emit:versions
script:
	def args1 = task.ext.args1?:""
	def args2 = task.ext.args2?:"--no-version --allow-overlaps --remove-duplicates "
	def args3 = optional_bed?"--regions-file \"${optional_bed}\"":""
	def limit = task.ext.limit?:10
	def prefix = task.ext.prefix?:(meta.id?:"variants")+".concat"
	def suffix = task.ext.suffix?:".bcf"
	def suffix2 = (suffix.endsWith("bcf")?"bcf":"vcf.gz")
	def suffix3 = (suffix.endsWith("bcf")?"bcf.csi":"vcf.gz.tbi")
	def by_contig = (task.ext.by_chromosome?:false) as boolean
"""	
	hostname 1>&2
	mkdir -p TMP
	find VCFS/ -name "*.vcf.gz" -o -name "*.bcf" > TMP/jeter.list
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
							${suffix.endsWith("bcf")?"--write-index":""} \\
							--threads ${task.cpus} \\
							${args1} \\
							${args2} \\
							--regions-file "TMP/\${CONTIG}.ctg.bed" \\
							-O ${suffix.endsWith("bcf")?"b9":"z9"} \\
							--file-list TMP/jeter2.vcf.list \\
							-o "TMP/jeter.${suffix2}" 

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
							${suffix.endsWith("bcf")?"--write-index":""} \\
							--threads ${task.cpus} \\
							${args1} \\
							${args2} \\
							--regions-file "TMP/\${CONTIG}.ctg.bed" \\
							-O ${suffix.endsWith("bcf")?"b9":"z9"} \\
							--file-list TMP/jeter3.list \\
							-o "TMP/jeter.${suffix2}" 

						rm -vf "TMP/*.delete.me.bcf"
						rm -vf "TMP/*.delete.me.bcf.csi"
					fi

					# default write index is CSI, not TBI for vcf.gz
					if ${!suffix.endsWith("bcf")}
					then
						bcftools index -f -t --threads ${task.cpus} TMP/jeter.${suffix2}
					fi

					mv -v TMP/jeter.${suffix2} "${prefix}.\${CONTIG}.${suffix2}"
					mv -v TMP/jeter.${suffix3} "${prefix}.\${CONTIG}.${suffix3}"

		done
		# end of loop over each chromosome


	else
		# do NOT group by chromosome
		
		if test  `wc -l < TMP/jeter.list` -le ${limit}
		then

			bcftools concat \\
				${suffix.endsWith("bcf")?"--write-index":""} \\
				--threads ${task.cpus} \\
				${args1} \\
				${args2} \\
				${args3} \\
				-O ${suffix.endsWith("bcf")?"b9":"z9"} \\
				--file-list TMP/jeter.list \\
				-o "TMP/jeter.${suffix2}" 

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
				${suffix.endsWith("bcf")?"--write-index":""} \\
				--threads ${task.cpus} \\
				${args1} \\
				${args2} \\
				${args3} \\
				-O ${suffix.endsWith("bcf")?"b9":"z9"} \\
				--file-list TMP/jeter2.list \\
				-o "TMP/jeter.${suffix2}" 

		fi

		# default write index is CSI, not TBI for vcf.gz
		if ${!suffix.endsWith("bcf")}
		then
			bcftools index -f -t --threads ${task.cpus} TMP/jeter.${suffix2}
		fi

		mv TMP/jeter.${suffix2} ${prefix}.${suffix2}
		mv TMP/jeter.${suffix3} ${prefix}.${suffix3}
	fi

cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
}
