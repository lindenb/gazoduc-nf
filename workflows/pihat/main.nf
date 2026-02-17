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


include {dumpParams;runOnComplete            } from '../../modules/utils/functions.nf'
include { BCFTOOLS_MERGE                           } from '../../modules/bcftools/merge3'
include {PIHAT as PIHAT01                          } from '../../subworkflows/pihat'
include {PREPARE_ONE_REFERENCE                     } from '../../subworkflows/samtools/prepare.one.ref'
include { MULTIQC                                  } from '../../subworkflows/multiqc'
include { VCF_INPUT                                } from '../../subworkflows/nf/vcf_input'
include { VCF_INPUT as VCF_INPUT_1KG               } from '../../subworkflows/nf/vcf_input'
include { VCF_INPUT as VCF_INPUT_GNOMAD            } from '../../subworkflows/nf/vcf_input'
include { META_TO_BAMS                             } from '../../subworkflows/samtools/meta2bams2'
include { READ_SAMPLESHEET                         } from '../../subworkflows/nf/read_samplesheet'
include { COMPILE_VERSIONS                         } from '../../modules/versions/main.nf'

if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}

workflow {
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
	if(params.vcf==null && params.samplesheet==null) {
		log.error("undefined option --vcf or --samplesheet")
		exit -1
		}
	else if(params.vcf!=null && params.samplesheet!=null) {
		log.error("both defined: --vcf or --samplesheet")
		exit -1
		}
	else if(params.vcf!=null) {
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
	    	.map{meta,vcf,tbi->[metadata,vcf,tbi]}
		}
	else 
		{
		if(params.gnomad==null) {
			log.error("--gnomad must be defined when input is bam")
			exit -1
			}

		/* Read samplesheet */
		READ_SAMPLESHEET(
			metadata.plus([arg_name:"samplesheet"]),
			params.samplesheet
			)
		versions = versions.mix(READ_SAMPLESHEET.out.versions)

		/** extract BAM/bai/fasta/fai/dict from samplesheet */
		META_TO_BAMS(
			metadata,
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
			PREPARE_ONE_REFERENCE.out.dict,
			READ_SAMPLESHEET.out.samplesheet
			)
		versions = versions.mix(META_TO_BAMS.out.versions)

		/** extract all frequent sites in gnomad  for each contig */
		FIND_SITES_IN_GNOMAD(
			PREPARE_ONE_REFERENCE.out.dict,
			PREPARE_ONE_REFERENCE.out.fai.splitCsv(sep:'\t',header:false)
				.map{meta,row->[meta,row[0]/* contig*/]}
				.filter{_meta,contig->contig.matches("(chr)?[0-9XY]+")}
				.combine(gnomad_ch)
				.map{_meta1,ctg,meta2,vcf,tbi->[meta2.plus(contig:ctg),vcf,tbi]}
			)
		versions = versions.mix(FIND_SITES_IN_GNOMAD.out.versions)

		/** genotype each BAM for each contig */
		GENOTYPE_GATK(
			PREPARE_ONE_REFERENCE.out.dict,
			META_TO_BAMS.out.bams.combine(FIND_SITES_IN_GNOMAD.out.vcf)
				.map{meta1,bam,bai,fasta,fai,dict,meta2,vcf,tbi->[meta1.plus(contig:meta2.contig),bam,bai,fasta,fai,dict,vcf,tbi]}
			)
		versions = versions.mix(GENOTYPE_GATK.out.versions)

		/** merge per contig */
		BCFTOOLS_MERGE(
			GENOTYPE_GATK.out.vcf
				.map{meta,vcf,tbi->[meta.contig,[vcf,tbi]]}
				.groupTuple()
				.map{contig,files->[ [id:contig],files.flatten().sort()]}
			)
		versions = versions.mix(BCFTOOLS_MERGE.out.versions)

		uservcf_ch = BCFTOOLS_MERGE.out.vcf
		}
	


	
	if(params.sample2pop==null) {
		sample2pop_ch = Channel.of([[id:"nosample2pop"],[]])
		}
	else
		{
		sample2pop_ch = Channel.of([[id:"sample2pop"],file(params.sample2pop)])
		}

	if(params.remove_samples==null) {
		exclude_samples_ch = Channel.of([[id:"no_x_samples"],[]])
		}
	else
		{
		exclude_samples_ch = Channel.of([[id:"exclude_samples"],file(params.remove_samples)])
		}
	if(params.exclude_bed==null) {
		exclude_bed_ch = Channel.of([[id:"no_x_bed"],[]])
		}
	else
		{
		exclude_bed_ch = Channel.of([[id:"exclude_bed"],file(params.exclude_bed)])
		}

	 PIHAT01(
		metadata.plus(level:1),
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
	multiqc = multiqc.mix(PIHAT01.out.multiqc)
	
	
	MULTIQC(
		metadata.plus("id":"pihat"),
		sample2pop_ch,
		versions,
		[[id:"no_mqc_config"],[]],
		multiqc
		)
	}

runOnComplete(workflow);


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
	def delta = (task.ext.delta?:0.4) as double // floriane fait à AF=10% 
	def min_af = task.ext.min_af?:(0.5 - delta)
	def max_af = task.ext.max_af?:(0.5 + delta)
	def jvm = task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"
"""
mkdir -p TMP
bcftools index --stats "${gnomad}" |\\
	awk -F '\t' '{C1=\$1; gsub(/^chr/,"",C1); C2=sprintf("chr%s",C1); if(C1=="${contig}" || C2=="${contig}") printf("%s\t0\t%d\\n",\$1,\$2);}' > TMP/jeter.bed

test -s TMP/jeter.bed

bcftools view -m 2 -M 2 -G --regions-file TMP/jeter.bed --types snps --apply-filters 'PASS,.' -i 'INFO/${tag} >= ${min_af} &&  INFO/${tag}<=${max_af}' -O u "${gnomad}"  |\\
	bcftools annotate -x 'ID,FILTER,INFO,QUAL' -O v |\\
	awk -F '\t' '/^#/{print;PREV=-1;C="";next;} {if(length(\$4)!=1 || length(\$5)!=1) next; if(\$1!=C || int(\$2)-PREV>${distance} ) {print;PREV=int(\$2);C=\$1;} }' |\\
	java ${jvm} -jar \${HOME}/jvarkit.jar vcfsetdict -R ${dict}  -n SKIP  |\\
	bcftools sort --max-mem '${task.memory.giga}G' -T TMP/tmp -O z -o TMP/jeter.vcf.gz


bcftools  index  --threads ${task.cpus} -f -t TMP/jeter.vcf.gz

mv TMP/jeter.vcf.gz ${prefix}.vcf.gz
mv TMP/jeter.vcf.gz.tbi ${prefix}.vcf.gz.tbi

touch versions.yml
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
