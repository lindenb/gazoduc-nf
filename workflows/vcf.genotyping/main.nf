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
include {BCFTOOLS_MERGE             }  from '../../modules/bcftools/merge3'


workflow {
		if(params.fasta==null) {
				throw new IllegalArgumentException("undefined --fasta");
				}
		if(params.samplesheet==null) {
				throw new IllegalArgumentException("undefined --samplesheet");
				}


		vcf= [[id:"vcf"],file(params.vcf),file(params.vcf+".tbi")]

		bams = Channel.fromPath(params.samplesheet)
			.splitCsv(header:true,sep:',')
			.map{row->[[id:row.sample],file(row.bam),file(row.bai)]}

		fasta = [[id:"reference"],file(params.fasta)]
		fai = [[id:"reference"],file(params.fasta+".fai")]
		dict = [[id:"reference"],file(params.fasta.replaceAll(".fasta",".dict"))]

		Channel.of(params.fasta).view()
		Channel.of(dict).view()

		if(params.method.equals("gatk")) {
			GATK_GENOTYPING(
				[id:"genotyping"],
				fasta,
				fai,
				dict,
				bams,
				vcf
				)
			}
		else if(params.method.equals("bcftools")) {
			//BCFTOOLS_GENOTYPING(genomeId,bams,vcf)
			}
		else
			{
			exit(-1,"unknown method ${params.method}");
			}
	}


workflow GATK_GENOTYPING {
	take:
		meta
		fasta
		fai
		dict
		bams
		vcf
	main:
		GENOTYPE_GATK(
			fasta,
			fai,
			dict,
			bams,
			vcf
			)
		BCFTOOLS_MERGE(GENOTYPE_GATK.out.vcf
			.map{meta,vcf,tbi->[vcf,tbi]}
			.flatMap()
			.collect()
			.map{files->[[id:"genotyping"],files.sort()]}
			)
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
		ch2 = BCFTOOLS_MERGEx(ch1.output.map{T->[T[0],T[1]+"\t"+T[2]]}.groupTuple())
	}



process GENOTYPE_GATK {
label "process_single"
tag "${meta.id}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
array 100
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta ),path(bam),path(bai)
	tuple val(meta4 ),path(vcf),path(tbi)

output:
	tuple val(meta),path("*.vcf.gz"),path("*.tbi"),emit:vcf
script:
	def prefix = task.ext.meta?:"${meta.id}"

"""
mkdir -p TMP


gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" HaplotypeCaller \\
	-I "${bam}" \\
	-L ${vcf} \\
	-R "${fasta}" \
	--alleles ${vcf}\\
	--output-mode EMIT_ALL_CONFIDENT_SITES \\
	-O "TMP/jeter2.vcf.gz"


bcftools annotate -x 'QUAL,^INFO/DP,INFO/AC,INFO/AN,INFO/AF' -O z -o "TMP/jeter3.vcf.gz" TMP/jeter2.vcf.gz
bcftools index -ft TMP/jeter3.vcf.gz



mv TMP/jeter3.vcf.gz ./${prefix}.vcf.gz
mv TMP/jeter3.vcf.gz.tbi  ./${prefix}.vcf.gz.tbi
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
	tuple val(L)
output:
	path("${params.prefix?:""}${contig}.merged.bcf"),emit:vcf
	path("${params.prefix?:""}${contig}.merged.bcf.csi"),emit:index
script:
"""
hostname 1>&2
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


process BCFTOOLS_MERGEx {
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
	def min_file_split = 75;
"""
hostname 1>&2
${moduleLoad("bcftools")}

mkdir -p TMP

cat << EOF | sort -T TMP -t '\t' -k1,1 | cut -f 2 > TMP/jeter.list
${L.join("\n")}
EOF

SQRT=`awk 'END{X=NR;if(${min_file_split} > 0 && X <= ${min_file_split}){print(X);} else {z=sqrt(X); print (z==int(z)?z:int(z)+1);}}' TMP/jeter.list`
split -a 9 --additional-suffix=.list --lines=\${SQRT} TMP/jeter.list TMP/chunck.

find TMP -type f -name "chunck.*.list" | sort | while read F
do
	echo "\${F}" 1>&2
	bcftools merge --threads ${task.cpus} --file-list "\${F}" --missing-to-ref  -O b -o "\${F}.bcf"
	bcftools index --threads ${task.cpus} "\${F}.bcf"
	echo "\${F}.bcf" >> TMP/jeter2.list
done

bcftools merge --threads ${task.cpus} --file-list TMP/jeter2.list --missing-to-ref  -O u |\
	bcftools +fill-tags -O b  -o "${params.prefix?:""}${contig}.merged.bcf"  -- -t  AN,AC,AF,AC_Hom,AC_Het,AC_Hemi,NS
bcftools index --threads ${task.cpus} "${params.prefix?:""}${contig}.merged.bcf"

"""
}
