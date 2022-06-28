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
include {moduleLoad} from './../../modules/utils/functions.nf'
include {SAMTOOLS_SAMPLES02} from './../../modules/samtools/samtools.samples.02.nf'

workflow GATK4_HAPCALLER_GVCFS_01 {
        take:
                meta
                reference
                bams
                beds
        main:
                version_ch = Channel.empty()

		samples_bams_ch = SAMTOOLS_SAMPLES02(meta.plus(["with_header":true]), reference, bams)
		version_ch = version_ch.mix(samples_bams_ch.version)


		each_bed = beds.splitText().map{it.trim()}
		each_bam = samples_bams_ch.output.splitCsv(header: true, sep : '\t')

		each_bam.dump()		

		hc_input_ch = each_bed.combine(each_bam).map{T->[
			"bed":T[0],
			"sample":T[1].sample,
			"bam":T[1].bam,
			"new_sample":T[1].new_sample,
			"reference": reference,
			"dbsnp": (meta.dbsnp?:""),
			"extraHC": (meta.extraHC?:""),
			"pedigree": (meta.pedigree?:""),
			"mapq": (meta.mapq?:"-1"),
			"bam_reference": T[1].reference
			]}
		hapcaller_ch = HAPCALLER_GVCF(meta, hc_input_ch)
		version_ch = version_ch.mix(hapcaller_ch.version)

	emit:
		version = version_ch
}

process HAPCALLER_GVCF {
tag "bed:${file(row.bed).name} ref:${file(row.bam_reference).name} bam:${file(row.bam).name}"
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
        tuple val(row),path("hapcaller.vcf.gz"),emit:bedvcf
        path("hapcaller.vcf.gz.tbi"),emit:index
        path("version.xml"),emit:version
script:
	def bam = row.bam?:""
	def sample = row.sample?:""
	def new_sample = row.new_sample?:""
	def extraHC = row.extraHC?:""
	def pedigree = row.pedigree?:""
	def mapq = row.mapq?:"-1"
	def bams = row.bams?:""
	def reference = row.reference?:""
	def bam_reference = row.bam_reference?:row.reference
	def bed = row.bed?:""
	def dbsnp = row.dbsnp?:""
"""
hostname 1>&2
${moduleLoad("gatk4 bcftools jvarkit")}

   set -o pipefail
   set -x
   mkdir TMP 

   if [ "${reference}" != "${bam_reference}" ] ; then
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/bedrenamechr.jar -R "${bam_reference}" --column 1 --convert SKIP "${bed}" > TMP/fixed.bed
	if [ ! -s TMP/fixed.bed ] ; then
		tail -n 1 "\${REF}.fai" | awk -F '\t' '{printf("%s\t0\t1\\n",\$1);}' > TMP/fixed.bed
   	fi 
   else
	ln -s "${bed}" TMP/fixed.bed
   fi


   gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" HaplotypeCaller \
     -I "${bam}" \
     -ERC GVCF \
     --seconds-between-progress-updates 600 \
     ${!dbsnp.isEmpty() &&  reference.equals(bam_reference)?"--dbsnp ${dbsnp}":""}   \
     -L TMP/fixed.bed \
     -R "${bam_reference}" \
     --sample-ploidy `awk -F '\t' 'BEGIN{P=1;} {if(!(\$1=="chrY" || \$1=="Y")) P=2;} END{print P;}' TMP/fixed.bed` \
     ${pedigree.isEmpty()?"":" --pedigree '${pedigree}'"} \
     ${(mapq as Integer)<1?"":" --minimum-mapping-quality '${mapq}'"} \
     -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
     ${extraHC} \
     -O "TMP/jeter.g.vcf.gz"

   # reference is not the main reference
   if [ "${reference}" != "${bam_reference}" ] ; then

	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/vcfsetdict.jar -n SKIP -R "${reference}" TMP/jeter.g.vcf.gz > TMP/jeter2.vcf
	
	bcftools sort -T TMP  --max-mem "${task.memory.giga}G" -O z -o TMP/jeter.g.vcf.gz TMP/jeter2.vcf
	bcftools index --force --tbi TMP/jeter.g.vcf.gz
	rm TMP/jeter2.vcf

	# annotate dbsnp
	if [ ! -z "${dbsnp}" ] ; then
		# bcftools annotate is just too slow

		gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"  VariantAnnotator \
			-R "${reference}" \
			--dbsnp "${dbsnp}" \
			-V TMP/jeter.g.vcf.gz \
			-O TMP/jeter2.g.vcf.gz
		mv TMP/jeter2.g.vcf.gz TMP/jeter.g.vcf.gz
		mv TMP/jeter2.g.vcf.gz.tbi TMP/jeter.g.vcf.gz.tbi
	fi

   fi
 
   # must rename sample
   if [ "${sample}" != "${new_sample}" ] ; then

	gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"  RenameSampleInVcf \
		INPUT=TMP/jeter.g.vcf.gz \
		OUTPUT=TMP/jeter2.g.vcf.gz \
		NEW_SAMPLE=${new_sample}

	mv TMP/jeter2.g.vcf.gz TMP/jeter.g.vcf.gz
	mv TMP/jeter2.g.vcf.gz.tbi TMP/jeter.g.vcf.gz.tbi
   fi



   mv TMP/jeter.g.vcf.gz  hapcaller.vcf.gz
   mv TMP/jeter.g.vcf.gz.tbi  hapcaller.vcf.gz.tbi


##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">haplotype caller in GVCF mode</entry>
        <entry key="bam">${bam}</entry>
        <entry key="sample">${sample}</entry>
        <entry key="new.sample">${new_sample}</entry>
        <entry key="reference">${reference}</entry>
        <entry key="bam_reference">${bam_reference}</entry>
        <entry key="dbsnp">${dbsnp}</entry>
</properties>
EOF

"""
}

