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
import {LINUX_SPLIT} from '../../modules/utils/split.nf'
import {BCFTOOLS_CONCAT_01} from '../bcftools/bcftools.concat.01.nf'

workflow GATK_HAPCALLER_DIRECT_01 {
	take:
		meta
		reference
		references
		bams
		beds
		dbsnp
		pedigree
	main:
		version_ch = Channel.empty()
		split_bams_ch = LINUX_SPLIT(["method":"--lines=${meta.nbams?:"10"}","suffix":".list"],bams)
		version_ch = version_ch.mix(split_bams_ch.version)
		
		each_bam_list = split_bams_ch.output.splitCsv(header: false,sep:'\t',strip:true).map{T->T[0]}
		capture_bed_ch = beds.splitCsv(header: false,sep:',',strip:true)
		bed_sample_bam_ch = capture_bed_ch.combine(capture_bed_ch)

		callbed_ch = CALL_BED(meta,reference,pedigree,bed_sample_bam_ch)
		version_ch = version_ch.mix(callbed_ch.version.first())

		ped2vcf_ch  = callbed_ch.bedvcf.groupTuple()
		mergedbed_ch = MERGE_BED(meta, reference, ped2vcf_ch)
		version_ch = version_ch.mix(mergedbed_ch.version)

		COLLECT_TO_FILE(meta,mergedbed_ch.collect())

		BCFTOOLS_CONCAT_01(meta,,"")

	emit:
		output = merge_ch.out
		version = version_ch
	}

process CALL_BED {
tag "${file(bed).name} ${file(bams_reference).name} ${file(bams).name}"
cache 'lenient'
memory {task.attempt==1?'10G':'10G'}
errorStrategy 'retry'
maxRetries 5
cpus 1
afterScript 'rm -rf TMP'
input:
	val(meta)
	val(reference)
	val(pedigree)
	val(dbsnp)
	tuple val(bed),val(bams_reference),val(bams)
output:
        tuple bed,path("genotyped.bcf"),emit:bedvcf
        path("genotyped.bcf.csi"),emit:index
        path("version.xml"),emit:version
script:
	def extraHC = meta.extraHC?:""
	def mapq = meta.mapq?:"-1"
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
</properties>
EOF
"""
}

