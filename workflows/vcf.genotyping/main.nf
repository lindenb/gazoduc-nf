include {runOnComplete;moduleLoad} from '../../modules/utils/functions.nf'
include {SAMTOOLS_SAMPLES} from '../../subworkflows/samtools/samtools.samples.03.nf'
include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf'

workflow {
	VCF_GENOTYPING(params.genomeId, file(params.bams), file(params.vcf))
	}

runOnComplete(workflow)

workflow VCF_GENOTYPING {
	take:
		genomeId
		bams
		vcf
	main:
		if(params.method.equals("gatk")) {
			GATK_GENOTYPING(genomeId,bams,vcf)
			}
		else if(params.method.equals("bcftools")) {
			BCFTOOLS_GENOTYPING(genomeId,bams,vcf)
			}
		else
			{
			exit(-1,"unknown method ${params.method}");
			}
	}


workflow GATK_GENOTYPING {
	take:
		genomeId
		bams
		vcf
	main:
		sn_bam_ch = SAMTOOLS_SAMPLES([:],bams)
		rows = sn_bam_ch.rows.
			filter{it.genomeId.equals(genomeId)}.
			map{it.plus("vcf":vcf)}
		ch1 = GENOTYPE_GATK(genomeId,rows)
		ch2= MERGE_VCF(ch1.output.splitCsv(header:false,sep:'\t').groupTuple())
	}

workflow BCFTOOLS_GENOTYPING {
	take:
		genomeId
		bams
		vcf
	main:
		v2b_ch = VCF_TO_BED([:],vcf)
		each_contig = v2b_ch.chromosomes.splitText().map{it.trim()}

		sn_bam_ch = SAMTOOLS_SAMPLES([:],bams)
		rows = sn_bam_ch.rows.
			filter{it.genomeId.equals(genomeId)}.
			map{it.plus("vcf":vcf)}.
			combine(each_contig).
			map{T->T[0].plus(contig:T[1])}
			
		ch1 = GENOTYPE_BCFTOOLS(genomeId,rows)
		ch2 = BCFTOOLS_MERGE(ch1.output.map{T->[T[0],T[1]+"\t"+T[2]]}.groupTuple())
	}



process GENOTYPE_GATK {
tag "${row.sample}"
afterScript "rm -rf TMP"
memory "5g"
input:
	val(genomeId)
	val(row)
output:
	path("output.tsv"),emit:output
script:
	def mapq = 30
	def bam = row.bam?:""
	def sample = row.sample?:""
	def reference = params.genomes[genomeId].fasta
	def dbsnp = "--dbsnp "+params.genomes[genomeId].dbsnp
"""
hostname 1>&2
mkdir -p TMP
${moduleLoad("gatk/0.0.0 bcftools")}

bcftools index -s "${row.vcf}" | cut -f 1 | while read C
do
	bcftools view -O z -o TMP/jeter.vcf.gz "${row.vcf}" "\${C}"
	bcftools index -ft TMP/jeter.vcf.gz

   gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" HaplotypeCaller \\
     -I "${row.bam}" \\
     ${dbsnp} \\
     -L TMP/jeter.vcf.gz \\
     -R "${reference}" \
     --minimum-mapping-quality '${mapq}' \\
     --alleles TMP/jeter.vcf.gz \\
     --output-mode EMIT_ALL_CONFIDENT_SITES \\
     -O "TMP/jeter2.vcf.gz"


     bcftools annotate -x 'QUAL,^INFO/DP,INFO/AC,INFO/AN,INFO/AF' -O b -o "TMP/${sample}.\${C}.bcf" TMP/jeter2.vcf.gz
     bcftools index -ft "TMP/${sample}.\${C}.bcf"


     rm -v TMP/jeter2.vcf.gz TMP/jeter2.vcf.gz.tbi TMP/jeter.vcf.gz TMP/jeter.vcf.gz.tbi


    echo "\${C}\t\${PWD}/${sample}.\${C}.bcf"  >> output.tsv
done

mv -v TMP/${sample}.*.bcf ./
mv -v TMP/${sample}.*.bcf.csi ./
"""
}


process GENOTYPE_BCFTOOLS {
tag "${row.sample} ${row.contig}"
afterScript "rm -rf TMP"
input:
	val(genomeId)
	val(row)
output:
	tuple val("${row.contig}"),val("${row.sample}"),path("${row.sample}.${row.contig}.bcf"),emit:output
script:
	def mapq = params.mapq?:30
	def bam = row.bam?:""
	def sample = row.sample?:""
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
	def ploidy = genome.ensembl_name
"""
hostname 1>&2
mkdir -p TMP
${moduleLoad("bcftools htslib")}

	bcftools query --regions "${row.contig}" -f'%CHROM\t%POS\t%REF,%ALT\\n' "${row.vcf}" |\\
		bgzip -c > TMP/jeter.tsv.gz

	tabix --force -s1 -b2 -e2 TMP/jeter.tsv.gz


	bcftools mpileup --redo-BAQ -a 'FORMAT/AD' -a 'FORMAT/DP' -a 'FORMAT/SP' -a 'FORMAT/SCR' -a 'FORMAT/QS' -e 'FORMAT/NMBZ' \\
		--threads ${task.cpus} \\
		--fasta-ref "${reference}" \\
		--regions-file "TMP/jeter.tsv.gz" -q ${mapq} -O u  -O u -o TMP/jeter2.bcf '${bam}'

	bcftools call  --keep-alts \\
		-a 'INFO/PV4' -a 'FORMAT/GQ' -a 'FORMAT/GP' \\
		--ploidy ${ploidy} \\
		--threads ${task.cpus}	\\
		--targets-file "TMP/jeter.tsv.gz" \\
		--constrain alleles \\
		--multiallelic-caller --output-type b -o "TMP/jeter3.bcf" TMP/jeter2.bcf


     mv -v TMP/jeter3.bcf "${sample}.${row.contig}.bcf"
     bcftools index --threads ${task.cpus} -f "${sample}.${row.contig}.bcf"

"""
}



process MERGE_VCF {
tag "${contig} N=${L.size()}"
memory "10g"
input:
	tuple 	val(contig),val(L)
output:
	path("${params.prefix?:""}${contig}.merged.bcf"),emit:vcf
	path("${params.prefix?:""}${contig}.merged.bcf.csi"),emit:index
script:
"""
hostname 1>&2
${moduleLoad("gatk bcftools")}
mkdir -p TMP
cat << EOF > TMP/jeter.list
${L.join("\n")}
EOF

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" MergeVcfs \\
     I=TMP/jeter.list \\
     O=TMP/jeter.vcf.gz

bcftools view -O b -o "${params.prefix?:""}${contig}.merged.bcf" TMP/jeter.vcf.gz
bcftools index "${params.prefix?:""}${contig}.merged.bcf"

"""
}


process BCFTOOLS_MERGE {
afterScript "rm -rf TMP"
tag "${contig} N=${L.size()}"
memory "10g"
cpus 10
input:
	tuple 	val(contig),val(L)
output:
	path("${params.prefix?:""}${contig}.merged.bcf"),emit:vcf
	path("${params.prefix?:""}${contig}.merged.bcf.csi"),emit:index
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}

mkdir -p TMP

cat << EOF | sort -T TMP -t '\t' -k1,1 | cut -f 2 > TMP/jeter.list
${L.join("\n")}
EOF


bcftools merge --threads ${task.cpus} --file-list TMP/jeter.list --missing-to-ref  -O u |\
	bcftools +fill-tags -O b  -o "${params.prefix?:""}${contig}.merged.bcf"  -- -t AN,AC,AF
bcftools index --threads ${task.cpus} "${params.prefix?:""}${contig}.merged.bcf"

"""
}
