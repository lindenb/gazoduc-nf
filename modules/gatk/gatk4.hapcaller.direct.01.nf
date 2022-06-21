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


process GATK4_HAPCALLER_DIRECT {
tag "${file(row.bed).name} ${file(row.bams_reference).name} ${file(row.bams).name}"
cache 'lenient'
memory {task.attempt==1?'10G':'20G'}
errorStrategy 'retry'
maxRetries 5
cpus 1
afterScript 'rm -rf TMP'
input:
	val(meta)
	val(row)
output:
        tuple val(row),path("genotyped.bcf"),emit:bedvcf
        path("genotyped.bcf.csi"),emit:index
        path("version.xml"),emit:version
script:
	def extraHC = row.extraHC?:""
	def mapq = row.mapq?:"-1"
	def pedigree = row.pedigree?:""
	def bams = row.bams?:""
	def reference = row.reference?:""
	def bams_reference = row.bams_reference?:reference
	def bed = row.bed?:""
	def dbsnp = row.dbsnp?:""
"""
hostname 1>&2
module load getModule("gatk4 bcftools jvarkit")


# allele specific annotation are not supported in non-gvcf mode
gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" HaplotypeCaller \
	-I "${bams}" \
	-L "${bed}" \
	-G StandardHCAnnotation -G StandardAnnotation \
	${extraHC} \
	--seconds-between-progress-updates 60 \
	${pedigree.isEmpty()?"":" -A PossibleDeNovo --pedigree "+pedigree} \
	${(mapq as Integer) > -1 ?"--minimum-mapping-quality "+mapq : "" } \
	${dbsnp.isEmpty() || !bams_reference.equals(reference)?"":"--dbsnp "+dbsnp} \
	-R "${bams_reference}" \
	-O "TMP/jeter.vcf.gz"

rm TMP/jeter.vcf.gz.tbi

# reference is not the main reference
if [ "${bams_reference}" != "${reference}" ] ; then

	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/vcfsetdict.jar -n SKIP -R "${reference}" TMP/jeter.vcf.gz > TMP/jeter2.vcf
	bcftools view --threads ${task.cpus} -O z -o TMP/jeter.vcf.gz TMP/jeter2.vcf
	bcftools index TMP/jeter.vcf.gz
	rm TMP/jeter2.vcf

	# annotate dbsnp
	if [ ! -z "${dbsnp}" ] ; then
		bcftools annotate --threads ${task.cpus}  -a "${dbsnp}" -c ID -O b -o TMP/jeter2.bcf TMP/jeter.vcf.gz
		bcftools view --threads ${task.cpus} -O z -o TMP/jeter.vcf.gz TMP/jeter2.bcf
		rm TMP/jeter2.bcf
	fi

fi


bcftools sort -T TMP  --max-mem "${task.memory.giga}G" -O b -o genotyped.bcf TMP/jeter.vcf.gz

bcftools index genotyped.bcf

###########################################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">call samples with gatk haplotypecaller. No GVCF</entry>
	<entry key="bams">${bams}</entry>
	<entry key="bams-count">\$(wc -l < "${bams}")</entry>
	<entry key="bed">${bed}</entry>
	<entry key="bed-count">\$(wc -l < "${bed}")</entry>
	<entry key="dbsnp">${dbsnp}</entry>
	<entry key="pedigree">${pedigree}</entry>
	<entry key="mapq">${mapq}</entry>
	<entry key="extraHC">${extraHC}</entry>
</properties>
EOF
"""
}

