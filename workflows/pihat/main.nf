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
nextflow.enable.dsl=2

include { validateParameters                       } from 'plugin/nf-schema'
include { paramsHelp                               } from 'plugin/nf-schema'
include { paramsSummaryLog                         } from 'plugin/nf-schema'
include { samplesheetToList                        } from 'plugin/nf-schema'
include { runOnComplete                            } from '../../modules/utils/functions'
include { parseBoolean                             } from '../../modules/utils/functions'
include { PIHAT                                    } from '../../subworkflows/pihat'
include { PREPARE_ONE_REFERENCE                    } from '../../subworkflows/samtools/prepare.one.ref'
include { MULTIQC                                  } from '../../subworkflows/multiqc'
include { VCF_INPUT                                } from '../../subworkflows/nf/vcf_input'
include { VCF_INPUT as VCF_INPUT_1KG               } from '../../subworkflows/nf/vcf_input'
include { VCF_INPUT as VCF_INPUT_GNOMAD            } from '../../subworkflows/nf/vcf_input'
include { COMPILE_VERSIONS                         } from '../../modules/versions/main.nf'



workflow {

	if(!workflow.stubRun) {
		validateParameters()
		}

	if( params.help ) {
		log.info(paramsHelp())
		exit 0
		}  else {
		// Print summary of supplied parameters
		log.info paramsSummaryLog(workflow)
		}



	versions = Channel.empty()
	multiqc = Channel.empty()
	def metadata = [
		id: "pihat"
		]
	
	PREPARE_ONE_REFERENCE(
		metadata.plus([skip_scatter:true]),
		Channel.fromPath(params.fasta).map{[[id:it.baseName],it]}
		)
    versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)


	/* load gnomad */
	if(params.gnomad==null) {
		gnomad_ch = Channel.empty()
		log.warn("No --gnomad defined")
		}
	else {
		VCF_INPUT_GNOMAD(metadata.plus([
			path: params.gnomad,
			arg_name: "gnomad",
			require_index : true,
			required: true,
			unique : false
			]))
		versions = versions.mix(VCF_INPUT_GNOMAD.out.versions)
		gnomad_ch = VCF_INPUT_GNOMAD.out.vcf
		}
	
	/* load 1000G */
	if(params.onekg==null) {
		onekg_ch = Channel.empty()
		log.warn("No --onekg defined")
		}
	else
		{
		/***************************************************
		*
		* VCF input
		*
		*/
		VCF_INPUT_1KG(metadata.plus([
			path: params.onekg,
			arg_name: "onekg",
			require_index : true,
			required: true,
			unique : false
			]))
		versions = versions.mix(VCF_INPUT_1KG.out.versions)
		onekg_ch = VCF_INPUT_1KG.out.vcf
		}

	/* load user's VCF */
	if(params.vcf==null ) {
		log.error("undefined option --vcf")
		exit -1
		}
	else  {
		/***************************************************
		*
		* VCF input
		*
		*/
		VCF_INPUT(metadata.plus([
			path: params.vcf,
			arg_name: "vcf",
			require_index : true,
			required: true,
			unique : false
			]))
		versions = versions.mix(VCF_INPUT.out.versions)
		uservcf_ch = VCF_INPUT.out.vcf
	    	.map{meta,vcf,tbi->[meta,vcf,tbi]}
		}
	
	/***************************************************
	*
	* samples to population
	*
	*/
	if(params.sample2pop==null) {
		sample2pop_ch = Channel.of([[id:"nosample2pop"],[]])
		}
	else
		{
		sample2pop_ch = Channel.of([[id:"sample2pop"],file(params.sample2pop)])
		}
	/***************************************************
	*
	* samples to exclude
	*
	*/
	if(params.remove_samples==null) {
		exclude_samples_ch = Channel.of([[id:"no_x_samples"],[]])
		}
	else
		{
		exclude_samples_ch = Channel.of([[id:"exclude_samples"],file(params.remove_samples)])
		}
	/***************************************************
	*
	* BED of regions to exclude
	*
	*/
	if(params.exclude_bed==null) {
		exclude_bed_ch = Channel.of([[id:"no_x_bed"],[]])
		}
	else
		{
		exclude_bed_ch = Channel.of([[id:"exclude_bed"],file(params.exclude_bed)])
		}

	 PIHAT(
		metadata.plus([
			level:1,
			contigs_regex:"${params.contigs_regex}"
			]),
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		sample2pop_ch,
		exclude_samples_ch,
		exclude_bed_ch,
		onekg_ch,
		gnomad_ch,
		uservcf_ch
		)
	multiqc = multiqc.mix(PIHAT.out.multiqc)
	versions = versions.mix(PIHAT.out.versions)
	
	if(parseBoolean(params.with_multiqc)) {
		MULTIQC(
			metadata.plus("id":"pihat"),
			sample2pop_ch,
			versions,
			[[id:"no_mqc_config"],[]],
			multiqc
			)
		}
	}

runOnComplete(workflow);


/**
 *
 * FIND Sites in gnomad that will be genotyped 
 *
 */
process FIND_SITES_IN_GNOMAD {
tag "${meta.id} ${meta.contig}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(dict)
	tuple val(meta ),path(gnomad),path(tbi)
output:
	tuple val(meta ),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
	path("versions.yml"),emit:versions
script:
	def tag=task.ext.tag?:"AF_nfe"
	def distance = task.ext.distance?:100
	def contig = "${meta.contig}"
	def prefix = task.ext.prefix?:"${meta.id}.${contig}"
	def delta_AF = (task.ext.delta_AF?:0.4) as double // floriane fait à AF=0.01
	def min_af = task.ext.min_af?:(0.5 - delta_AF)
	def max_af = task.ext.max_af?:(0.5 + delta_AF)
	def jvm = task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"
"""
mkdir -p TMP
bcftools index --stats "${gnomad}" |\\
	awk -F '\t' '{C1=\$1; gsub(/^chr/,"",C1); C2=sprintf("chr%s",C1); if(C1=="${contig}" || C2=="${contig}") printf("%s\t0\t%d\\n",\$1,\$2);}' > TMP/jeter.bed

test -s TMP/jeter.bed

bcftools view -m 2 -M 2 -G --regions-file TMP/jeter.bed --types snps --apply-filters 'PASS,.' -i 'INFO/${tag} >= ${min_af} &&  INFO/${tag}<=${max_af} && INFO/NEGATIVE_TRAIN_SITE==0' -O u "${gnomad}"  |\\
	bcftools annotate -x 'ID,FILTER,INFO,QUAL' -O u |\\
	bcftools norm --multiallelics +any -O u |\\
	bcftools view -m 2 -M 2 -O v |\\
	awk -F '\t' '/^#/{print;PREV=-1;C="";next;} {if(length(\$4)!=1 || length(\$5)!=1) next; if(\$1!=C || int(\$2)-PREV>${distance} ) {print;PREV=int(\$2);C=\$1;} }' |\\
	java ${jvm} -jar \${HOME}/jvarkit.jar vcfsetdict -R ${dict}  -n SKIP  |\\
	bcftools sort --max-mem '${task.memory.giga}G' -T TMP/tmp -O z -o TMP/jeter.vcf.gz


bcftools  index  --threads ${task.cpus} -f -t TMP/jeter.vcf.gz

mv TMP/jeter.vcf.gz ${prefix}.vcf.gz
mv TMP/jeter.vcf.gz.tbi ${prefix}.vcf.gz.tbi

touch versions.yml
"""
stub:
	def contig = "${meta.contig}"
	def prefix = task.ext.prefix?:"${meta.id}.${contig}"
"""
touch versions.yml  ${prefix}.vcf.gz ${prefix}.vcf.gz.tbi
"""
}

process GENOTYPE_GATK {
label "process_single"
tag "${meta.id} ${meta.contig}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
        tuple val(meta1), path("_ref.dict")
        tuple val(meta ), path(bam),path(bai),path(fasta),path(fai),path(dict),path(vcf),path(tbi)
output:
        tuple val(meta /* contains contig */),path("*.vcf.gz"),path("*.tbi"),emit:vcf
		path("versions.yml"),emit:versions
script:
        def prefix = task.ext.meta?:"${meta.id}.${meta.contig}"
		def jvm = task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"
"""
mkdir -p TMP
set -x
if java ${jvm} -jar \${HOME}/jvarkit.jar samedict ${dict} ${vcf}
then
	cp -v "${vcf}" TMP/genotypes.vcf.gz
	cp -v "${tbi}" TMP/genotypes.vcf.gz.tbi
else
	java ${jvm} -jar \${HOME}/jvarkit.jar vcfsetdict -R "${fasta}"  -n SKIP ${vcf} |\\
	 	bcftools sort --max-mem '${task.memory.giga}G' -T TMP/tmp -O z -o TMP/genotypes.vcf.gz
	
	bcftools index  --threads ${task.cpus} -f -t TMP/genotypes.vcf.gz
fi

gatk --java-options "${jvm}" HaplotypeCaller \\
        -I "${bam}" \\
        -L  TMP/genotypes.vcf.gz \\
        -R "${fasta}" \\
        --alleles TMP/genotypes.vcf.gz\\
        --output-mode EMIT_ALL_CONFIDENT_SITES \\
        -O "TMP/jeter2.vcf.gz"


if java ${jvm} -jar \${HOME}/jvarkit.jar samedict _ref.dict TMP/jeter2.vcf.gz
then
	mv TMP/jeter2.vcf.gz TMP/jeter3.vcf.gz
else
	java ${jvm} -jar \${HOME}/jvarkit.jar vcfsetdict -R _ref.dict  -n SKIP TMP/jeter2.vcf.gz |\\
	 	bcftools sort --max-mem '${task.memory.giga}G' -T TMP/tmp -O z -o TMP/jeter3.vcf.gz
fi

bcftools annotate --threads  ${task.cpus}  -x 'ID,QUAL,FILTER,^INFO/DP,^FORMAT/GT' -O z -o TMP/jeter2.vcf.gz TMP/jeter3.vcf.gz
bcftools index --threads ${task.cpus} -f -t TMP/jeter2.vcf.gz

mv TMP/jeter2.vcf.gz ./${prefix}.vcf.gz
mv TMP/jeter2.vcf.gz.tbi  ./${prefix}.vcf.gz.tbi
touch versions.yml
"""

stub:
        def prefix = task.ext.meta?:"${meta.id}"
"""
	touch ${prefix}.vcf.gz ${prefix}.vcf.gz.tbi versions.yml
"""
}
