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
process GATK_BAM2VCF {
tag "${bed.name}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path("_reference.fa")
	tuple val(meta2),path("_reference.fa.fai")
	tuple val(meta3),path("_reference.dict")
	tuple val(meta4),path(dbsnp),path(dbsnp_tbi)
	tuple val(meta5),path(pedigree)
	tuple val(meta6),path("REFS/*")	
	tuple val(meta),path("BAMS/*"),path(bed)
output:
	tuple  val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),path(bed),emit:vcf
	path("versions.yml"),emit:versions
script:
	def fasta = "_reference.fa";
	def fai = "${fasta}.fai";
	def min_file_split = task.ext.min_file_split?:25
	def args1 = task.ext.args1?:""
	def args2 = task.ext.args2?:""
	def args3 = task.ext.args3?:""
	def args4 = task.ext.args4?:""
	def prefix = task.ext.prefix?:"${bed.name}"
	def jvm = "-Xmx${task.memory.giga}g -XX:-UsePerfData -Djava.io.tmpdir=TMP -DGATK_STACKTRACE_ON_USER_EXCEPTION=true"
"""
hostname 1>&2
mkdir -p TMP
set -x

touch TMP/samples.txt
echo "${fasta}" > TMP/references.list
find \${PWD}/REFS/ \\( -name "*.fa" -o  -name "*.fna" -o -name "*.fasta" \\) | sort -T TMP -V >> TMP/references.list

find ./BAMS/ \\( -name "*.bam" -o -name "*.cram" \\) | sort -T TMP -V > TMP/bams.list

i=1
cat TMP/bams.list | samtools samples -F  TMP/references.list | while read SAMPLE BAM FASTA
do

	# test fasta is known
	test "\${FASTA}" != "."

	# no the same reference ? change the BED according to chr notation
	if cmp "\${FASTA}.fai" "${fai}"
	then
		cp -v "${bed}" TMP/${bed.name}
	else
		jvarkit ${jvm} bedrenamechr \\
			-f "\${FASTA}" --column 1 --convert SKIP "${bed}" > TMP/${bed.name}
	
		if ! test -s TMP/${bed.name}
		then
			awk -F '\t' '{printf("%s\t0\t1\\n",\$1);}' "\${FASTA}.fai" | tail -n 1 > TMP/${bed.name}
		fi
	fi


   gatk --java-options "${jvm}" HaplotypeCaller \\
     -L "TMP/${bed.name}" \\
     -R "\${FASTA}" \\
     -I "\${BAM}" \\
     -ERC GVCF \\
     ${args1} \\
     -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \\
     -O "TMP/jeter.\${i}.g.vcf.gz"

     rm -f TMP/tmp_read_resource_*.config TMP/libgkl_*.so

	
	if ! cmp "\${FASTA}.fai" "${fai}"
	then
			jvarkit ${jvm} vcfsetdict \\
					-n SKIP \\
					-R "${fasta}" \\
					"TMP/jeter.\${i}.g.vcf.gz" > TMP/jeter2.vcf

			bcftools sort  -T TMP/sort  --max-mem "${task.memory.giga}G" -O z -o "TMP/jeter.\${i}.g.vcf.gz" TMP/jeter2.vcf
			bcftools index --threads ${task.cpus} --force --tbi "TMP/jeter.\${i}.g.vcf.gz"
	fi


	SAMPLE2="\${SAMPLE}"
	while grep -Fxq "\${SAMPLE2}" TMP/samples.txt
	do
		SAMPLE2="\${SAMPLE2}.x"
	done

	if test "\${SAMPLE2}" != "\${SAMPLE}"
	then
		gatk --java-options "${jvm}"  RenameSampleInVcf \\
				-INPUT "TMP/jeter.\${i}.g.vcf.gz" \\
				-OUTPUT TMP/jeter2.g.vcf.gz \\
				-NEW_SAMPLE_NAME "\${SAMPLE2}"

		bcftools index  --threads ${task.cpus}  --force --tbi TMP/jeter2.g.vcf.gz
		mv TMP/jeter2.g.vcf.gz "TMP/jeter.\${i}.g.vcf.gz"
		mv TMP/jeter2.g.vcf.gz.tbi "TMP/jeter.\${i}.g.vcf.gz"
	fi



	echo "\${SAMPLE2}" >>  TMP/samples.txt

    echo "TMP/jeter.\${i}.g.vcf.gz" >> TMP/gvcfs.list
    i=\$((i+1))		

done

SQRT=`awk 'END{X=NR;if(X <= ${min_file_split}){print(X);} else {z=sqrt(X); print (z==int(z)?z:int(z)+1);}}' TMP/gvcfs.list`

split -a 9 --additional-suffix=.list --lines=\${SQRT} TMP/gvcfs.list  TMP/split.
rm TMP/gvcfs.list

i=1
find ./TMP -type f -name "split*.list" | while read F
do
	gatk --java-options "${jvm}" \\
		CombineGVCFs \\
		-R "${fasta}" \\
		-L "${bed}" \\
		-V "\${F}" \\
		${args2} \\
		-O "TMP/combined.\${i}.g.vcf.gz" \\
		-G StandardAnnotation \\
		-G AS_StandardAnnotation

	 echo "TMP/combined.\${i}.g.vcf.gz" >> TMP/gvcfs.list
	 i=\$((i+1))
  
         rm -f TMP/tmp_read_resource_*.config TMP/libgkl_*.so

done

if [[ \$(wc -l < TMP/gvcfs.list) -gt 1  ]]
then

	gatk --java-options "${jvm}" \\
		CombineGVCFs \\
		-R "${fasta}" \\
		-L "${bed}" \\
		${args3} \\
		-V "TMP/gvcfs.list" \\
		-O TMP/combined.all.g.vcf.gz \\
		-G StandardAnnotation \\
		-G AS_StandardAnnotation

else

	mv -v TMP/combined.1.g.vcf.gz     TMP/combined.all.g.vcf.gz
	mv -v TMP/combined.1.g.vcf.gz.tbi TMP/combined.all.g.vcf.gz.tbi

fi



gatk --java-options "${jvm}" \\
	GenotypeGVCFs \\
	-R "${fasta}" \\
	-L "${bed}" \\
        ${dbsnp?"--dbsnp \"${dbsnp}\"":""} \\
        ${pedigree?"--pedigree \"${pedigree}\"":""} \\
	${args4} \\
	-O TMP/genotyped.vcf.gz \\
	-V TMP/combined.all.g.vcf.gz \\
	-G StandardAnnotation \\
	-G AS_StandardAnnotation



mv -v TMP/genotyped.vcf.gz "${prefix}.vcf.gz"
mv -v TMP/genotyped.vcf.gz.tbi "${prefix}.vcf.gz.tbi"

cat << EOF > versions.yml
${task.process}:
    gatk: "\$( gatk --version 2> /dev/null  | paste -s -d ' ' )"
EOF
"""

stub:
def prefix = task.ext.prefix?:"${bed.name}"
"""
touch versions.yml "${prefix}.vcf.gz" "${prefix}.vcf.gz.tbi"
"""
}
