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
nextflow.enable.dsl=2

include {k1_signature} from '../../modules/utils/k1.nf'
include {runOnComplete;dumpParams} from '../../modules/utils/functions.nf'
include {MAKE_STATS as MAKE_STATS_MALE;MAKE_STATS as MAKE_STATS_FEMALE; MAKE_STATS as MAKE_STATS_BOTH} from  './sub.nf'
if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}


workflow {
	par_ch = PAR_TO_BED(file(params.fai))
	xx_fasta_ch = MAKE_FEMALE_FASTA(file(params.fasta), file(params.fai), file(params.dict) )

	sex_fa_fai_dict_ch = xx_fasta_ch.map{["XX",it[0],it[1],it[2]]}.
		mix(Channel.of([file(params.fasta),file(params.fai),file(params.dict)]).map{["XY",it[0],it[1],it[2]]})

	bwa_idx_ch = BWA_INDEX_GENOME(sex_fa_fai_dict_ch)
	
	ch1 = Channel.fromPath(params.samplesheet).
		splitCsv(header:true, sep:'\t').
		map{T->{
			if(T.sex.equals("female")) return T.plus("sex":"XX");
			if(T.sex.equals("male")) return T.plus("sex":"XY");
			}}.
		map{[it.sample,it.bam,it.bai,it.sex,it.fasta,it.fai]}.
		combine(bwa_idx_ch).
		filter{it[3].equals(it[6])}
	
	remap_ch = REMAP_XY(ch1)

	contigs_ch = Channel.fromPath(params.fai).
		splitCsv(header:false,sep:'\t').
		map{it[0]}.
		filter{it.matches("(chr)?[XY]")}

	par2_ch = par_ch.par.map{["par",it]}.mix(par_ch.non_par.map{["hap",it]})

	par3_ch = remap_ch.output.
		combine(par2_ch).
		combine(contigs_ch)
		
		// NO we need the GVCF stats filter{!(it[3].equals("XX") && it[9].matches("(chr)?Y"))}


	gvcf_ch = GATK_GVCF(par3_ch)
	ch3 = gvcf_ch.output.
		map{[[it[0],it[1]],it[2]]}.
		groupTuple().
		map{[it[0],it[1].flatten()]}

	windows_ch = MAKE_WINDOWS(par2_ch)
	intervals_ch = windows_ch.output.splitCsv(header:false,sep:'\t').
		map{[it[0][0],it[0][0]+":"+(1+(it[0][1] as int))+"-"+it[0][2],it[1]]}.
		combine(ch3).
		filter{it[0].equals(it[3][0]) /* chromosome */ && it[2].equals(it[3][1])/* hap */}.
		map{[it[1] /* intervals */, it[4] /* vcfs */]}

	
	ch4 = GATK_GENOTYPE(
		file(params.fasta),
		file(params.fai),
		file(params.dict),
		intervals_ch
		)
	
	ch5 = CONCAT_VCFS(ch4.flatten().collect())
	
	MAKE_STATS_MALE(
		file(params.fasta),
		file(params.fai),
		file(params.gtf),
		file(params.samplesheet),
		ch5.output,
		"M"
		)
	MAKE_STATS_FEMALE(
		file(params.fasta),
		file(params.fai),
		file(params.gtf),
		file(params.samplesheet),
		ch5.output,
		"F"
		)
	MAKE_STATS_BOTH(
		file(params.fasta),
		file(params.fai),
		file(params.gtf),
		file(params.samplesheet),
		ch5.output,
		"MF"
		)
	}
	

runOnComplete(workflow)

process PAR_TO_BED {
tag "${fai.name}"
label "process_quick"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	path(fai)
output:
	path("par.bed"),emit:par
	path("non.par.bed"),emit:non_par
script:
	def k1 = k1_signature()
"""
mkdir -p TMP
set -o pipefail

export LC_ALL=C

cat << __EOF__  > TMP/jeter.csv
1:${k1.hg38},X:10001-2781479
1:${k1.hg38},Y:10001-2781479
1:${k1.hg38},X:155701383-156030895
1:${k1.hg38},Y:56887903-57217415
1:${k1.hg19},X:60001-2699520
1:${k1.hg19},Y:10001-2649520
1:${k1.hg19},X:154931044-155260560
1:${k1.hg19},Y:59034050-59363566
__EOF__

# duplicate with 'chr' suffix
awk -F, '{printf("%s\\nchr%s,chr%s\\n",\$0,\$1,\$2);}' TMP/jeter.csv |\\
	sort -T TMP -t, -k1,1 > TMP/jeter2.csv

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' |\\
	sort -T TMP -t, -k1,1 > TMP/jeter3.csv

join -t, -1 1 -2 1 -o '2.2' TMP/jeter3.csv TMP/jeter2.csv |\\
	tr ":-" "\\t" |\\
	sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/par.bed

test -s TMP/par.bed

cut -f1,2 '${fai}' |\\
	sort -T TMP -t, -k1,1 -k2,2n > TMP/genomed.bed

bedtools complement -i TMP/par.bed -L -g TMP/genomed.bed |\\
	sort -T TMP -t, -k1,1 -k2,2n > TMP/non.par.bed

test -s TMP/non.par.bed

mv TMP/non.par.bed ./
mv TMP/par.bed ./
"""
}

process MAKE_FEMALE_FASTA {
tag "${fasta.name}"
label "process_quick"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	path(fasta)
	path(fai)
	path(dict)
output:
	tuple path("female.${fasta.simpleName}.fa"), path("female.${fasta.simpleName}.fa.fai"), path("female.${fasta.simpleName}.dict"),emit:output
script:
"""
mkdir -p TMP
set -o pipefail


awk -F '\t' '\$1 ~ /^(chr)?Y\$/ {printf("%s\t0\t%s\\n",\$1,\$2);}' "${fai}" > TMP/mask.bed
test -s TMP/mask.bed

bedtools maskfasta -fi "${fasta}" -bed TMP/mask.bed -fo TMP/jeter.fa

mv TMP/jeter.fa female.${fasta.simpleName}.fa

## NO because is the original fasta file, there is can be comments after the name of the chromosome
## cp ${fai} female.${fasta.simpleName}.fa.fai

samtools faidx -o female.${fasta.simpleName}.fa.fai female.${fasta.simpleName}.fa 

cp ${dict} female.${fasta.simpleName}.dict
"""
}

process BWA_INDEX_GENOME {
tag "${fasta.name}"
label "process_quick"
conda "${moduleDir}/../../conda/bwa.yml"
afterScript "rm -rf TMP"
input:
	tuple val(name),path(fasta),path(fai),path(dict)
output:
    tuple val(name),path("${name}"),path(fasta),path(fai),path(dict),emit:output
script:
"""
    mkdir -p TMP
    bwa index -p TMP/${name} ${fasta}
    mv TMP "${name}"
"""
}

process REMAP_XY {
tag "${sample} ${bam.name} ${sex}"
label "process_quick_high"
conda "${moduleDir}/../../conda/bwa.yml"
afterScript "rm -rf TMP"

input:
	tuple val(sample),path(bam),path(bai),val(sex),path("bamref.fasta"),path("bamref.fasta.fai"),val(ignore),path(bwa_dir),path(fasta),path(fai),path(dict)
output:
	tuple val(sample),path("${sample}.bam"), path("${sample}.bam.bai"),val(sex),path(fasta),path(fai),path(dict),emit:output
script:
"""
mkdir -p TMP
set -o pipefail
set -x
BWA_INDEX=`find ${bwa_dir}/ -name "*.amb" | sed 's/\\.amb\$//'`
test ! -z "\${BWA_INDEX}"

samtools view --fetch-pairs --reference "bamref.fasta" --uncompressed -O BAM "${bam}" `awk -F '\t' '\$1 ~ /^(chr)?(X|Y)\$/' "${fai}" | cut -f1 | paste -sd ' '` |\\
	samtools collate -f -O -u --no-PG --reference "bamref.fasta" - TMP/tmp.collate |\\
	samtools fastq -N -1 TMP/jeter.R1.fq.gz -2 TMP/jeter.R2.fq.gz -s TMP/jeter.Rx.fq.gz -0 TMP/jeter.R0.fq.gz -n

bwa mem -t ${task.cpus} \\
		-R '@RG\\tID:${sample}\\tSM:${sample}\\tLB:${sample}\\tCN:BirdNantes\\tPL:ILLUMINA' \\
		"\${BWA_INDEX}" TMP/jeter.R1.fq.gz  TMP/jeter.R2.fq.gz |\\
		samtools view --uncompressed -O BAM -o TMP/jeter.bam '-'

	# collate
	samtools collate -l 5 \\
		--threads ${task.cpus} \\
		--output-fmt BAM \\
		--no-PG  \\
		-T TMP/tmp2.collate \\
		-o TMP/jeter2.bam \\
		TMP/jeter.bam
	mv TMP/jeter2.bam TMP/jeter.bam

	# fixmate
	samtools fixmate \\
		--threads  ${task.cpus} \\
		-mc \\
		--output-fmt BAM \\
		TMP/jeter.bam \\
		TMP/jeter2.bam
	mv TMP/jeter2.bam TMP/jeter.bam

	# sort
	samtools sort \\
		--threads ${task.cpus} \\
		-o TMP/jeter2.bam \\
		-O BAM \\
		-T TMP/tmp \\
		TMP/jeter.bam
	mv TMP/jeter2.bam TMP/jeter.bam

	# markdup
	samtools markdup \\
		-s --json -f TMP/stats.json \\
		-T TMP/tmp \\
		--threads ${task.cpus} \\
		TMP/jeter.bam TMP/jeter2.bam
	mv TMP/jeter2.bam TMP/jeter.bam

	# index
	samtools index --threads ${task.cpus}  TMP/jeter.bam

mv TMP/jeter.bam ${sample}.bam
mv TMP/jeter.bam.bai ${sample}.bam.bai
mv TMP/stats.json ${sample}.stats.json
"""
}


process GATK_GVCF {
tag "${sample} ${par_bed.name} ${contig} ${par} ${fasta.name} ${sex}"
label "process_quick"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
array 100
afterScript "rm -rf TMP"
input:
	tuple val(sample), path(bam), path(bai), val(sex),path(fasta), path(fai), path(dict), val(par), path(par_bed), val(contig)
output:
	tuple val(contig),val(par),path("${sample}.${contig}.*"),emit:output
script:
	def ploidy = ((sex.equals("XX") && contig.matches("(chr)?X")) || par.equals("par")?2:1)
"""
mkdir -p TMP

   awk -F '\t' '(\$1=="${contig}")' "${par_bed}"  > TMP/select.bed

   gatk --java-options "-XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" HaplotypeCaller \\
     -I "${bam}" \\
     -ERC GVCF \\
     --sample-ploidy ${ploidy} \\
     -L TMP/select.bed \\
     -R "${fasta}" \\
     -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \\
     -O "TMP/jeter.g.vcf.gz"

mv TMP/jeter.g.vcf.gz ${sample}.${contig}.${par}.vcf.gz
mv TMP/jeter.g.vcf.gz.tbi ${sample}.${contig}.${par}.vcf.gz.tbi
"""
}


process MAKE_WINDOWS {
tag "${type}"
label "process_quick"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(type),path(bed)
output:
	tuple path("*.bed"),val(type),emit:output
script:
	def w=1000000
	def s=w-1000
"""
mkdir -p TMP
set -o pipefail

sort -T TMP -t '\t' -k1,1 -k2,2n "${bed}" |\\
	bedtools makewindows -b - -w ${w} -s ${s} > TMP/jeter.bed

mv TMP/jeter.bed ${type}_${w}_${s}.bed
"""
}

process GATK_GENOTYPE {
tag "${interval}"
label "process_quick"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
        path(fasta)
	path(fai)
	path(dict)
        tuple val(interval),path("VCFS/*")
output:
        path("genotyped.*"),emit:output
script:
	def title = "genotyped." + interval.replace(':','_').replace('-','_')
"""
hostname 1>&2
mkdir -p TMP
set -x

find ./VCFS -name "*.vcf.gz" > TMP/jeter.list
test -s TMP/jeter.list

cat TMP/jeter.list | while read F
do
        gunzip -c "\${F}" | grep "^#CHROM" -m1 | cut -f 10 | tr "\\n" "\t" >> TMP/sample.map
        echo "\${F}" >> TMP/sample.map
done


# sort on sample name
LC_ALL=C sort -t '\t' -k1,1 -T TMP TMP/sample.map > "TMP/jeter.map"
test -s TMP/jeter.map
mv TMP/jeter.map TMP/sample.map

SQRT=`awk 'END{X=NR;if(X<10){print(X);} else {z=sqrt(X); print (z==int(z)?z:int(z)+1);}}' TMP/sample.map`

gatk --java-options "-XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" GenomicsDBImport \\
    -R ${fasta} \\
    --batch-size \${SQRT} \\
    ${(task.cpus as int) > 1 ? "  --reader-threads " +task.cpus : "" } \\
    --sample-name-map TMP/sample.map \\
    -L "${interval}" \\
    --genomicsdb-workspace-path "TMP/database"

gatk --java-options "-XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" GenotypeGVCFs  \\
     -R "${fasta}" \\
      -L "${interval}" \\
      -V gendb://TMP/database \\
      -O "TMP/jeter.vcf.gz"

mv TMP/jeter.vcf.gz "${title}.vcf.gz"
mv TMP/jeter.vcf.gz.tbi "${title}.vcf.gz.tbi"
"""
}

process CONCAT_VCFS {
label "process_quick"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	path("VCFS/*")
output:
	path("gatk4.xy.*"),emit:output
script:
"""
hostname 1>&2
mkdir -p TMP

find ./VCFS -name "*.vcf.gz" > TMP/jeter.list
test -s TMP/jeter.list

bcftools concat -a --remove-duplicates --file-list TMP/jeter.list -O b9 --threads ${task.cpus} -o TMP/gatk4.xy.bcf 
bcftools index -f TMP/gatk4.xy.bcf

mv TMP/gatk4.xy.* ./
"""
}
