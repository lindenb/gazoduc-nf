/*

Copyright (c) 2024 Pierre Lindenbaum

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

include {BWA_INDEX} from '../../../modules/bwa/index'


boolean isEmptyGz(Object filename) {
	final java.nio.file.Path p;
	if(filename instanceof java.nio.file.Path) {
		p = (java.nio.file.Path)filename;
		}
	else
		{
		p = java.nio.file.Paths.get(filename.toString());
		}
	if(!java.nio.file.Files.exists(p)) return true;
	try (java.util.zip.GZIPInputStream in = new  java.util.zip.GZIPInputStream(java.nio.file.Files.newInputStream(p)) ) {
		return in.read()==-1;
		}
	}

workflow REMAP_BWA_01 {
	take:
		genome
		samplesheet
		bed
	main:
		version_ch = Channel.empty()
	
		bams_ch = samplesheet.splitCsv(header:true,sep:'\t').
			map{T->{
				if(T.sample==null || T.sample.isEmpty()) throw new IllegalArgumentException("sample missing");
				def sample = T.sample;
				if(T.bam==null || T.bam.isEmpty()) throw new IllegalArgumentException("bam missing");
				def bam = T.bam;
				def bai = (T.bai == null ?(bam.endsWith(".bam")?bam+".bai":bam+".crai") : T.bai);
				def fasta = "NO_FASTA_IN";
				def fai = "NO_FAI_IN";
				if(T.bam.endsWith(".cram")) {
					if(T.fasta==null || T.fasta.isEmpty()) throw new IllegalArgumentException("fasta missing for CRAM");
					fasta  = T.fasta
					fai = (T.fai==null || T.fai.isEmpty() ? fasta+".fai": T.fai)
					}
				def bed = (T.bed==null || T.bed.isEmpty() || T.bed.equals(".")?"NO_BED":T.bed);
				return [sample,file(bam),file(bai),file(fasta),file(fai),file(bed)];
				}}

		index_bwa = BWA_INDEX(file(params.fasta))
		remap_ch = REMAP_ONE(genome, bed, index_bwa.output ,bams_ch )

	emit:
		output = index_bwa.output
	}

workflow REMAP_ONE {
take:
	genome
	bed
	index_bwa
	bams_ch
main:
	if(params.collate==true) {
	        collate_ch = SAMTOOLS_COLLATE(bams_ch)
		ch1 = collate_ch.output.flatMap{T->
			[
			[T[0],T[1],T[2]],
			[T[0],T[3],file("NO_FILE")],
			[T[0],T[4],file("NO_FILE")]
			]}.
			filter{!(isEmptyGz(it[1]) && isEmptyGz(it[2]))}
		map_ch = BWA_MEM_01(genome, index_bwa, bed, ch1)
		}
	else
		{
		map_ch = COLLATE_AND_MAP_BWA(genome,index_bwa,bed,bams_ch)
		}
	mapped_ch = MERGE_BAMS(genome, map_ch.output.groupTuple().map{[it[0],it[1].flatten().collect()]})
	

emit:
	bams = mapped_ch.output
}



process SAMTOOLS_COLLATE {
tag "${sample} ${bam.name}"
label "process_medium"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(sample),path(bam),path(bai),path(fasta),path(fai),path(bed)
output:
	tuple val(sample),
			path("FASTQS/${sample}.paired.R1.fq.gz"),
			path("FASTQS/${sample}.paired.R2.fq.gz"),
			path("FASTQS/${sample}.singleton.fq.gz"),
			path("FASTQS/${sample}.other.fq.gz"),emit:output
script:
	def fetch_pairs = task.attempt==1?" --fetch-pairs ":""
	def level = task.ext.compression_level?:1
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP


# extract reads
if ${bed.name.contains(".")} ; then
	samtools view  \\
		--fast \\
		--threads ${task.cpus} ${fetch_pairs} \\
		--reference "${fasta}" \\
		-M -O BAM -o TMP/jeter.bam \\
		--regions-file "${bed}" \\
		"${bam}"
fi


# collate
samtools collate -f \\
	--threads ${(task.cpus as int) -1} \\
	-O -u --no-PG \\
	--reference "${fasta}" \\
	"${bed.name.contains(".") ? "TMP/jeter.bam":"${bam}"}" TMP/tmp.collate |\
	samtools fastq -n --threads 1 \\
		-1 TMP/${sample}.paired.R1.fq.gz \\
		-2 TMP/${sample}.paired.R2.fq.gz \\
		-s "TMP/${sample}.singleton.fq.gz" \\
		-0 "TMP/${sample}.other.fq.gz"


if ${bed.name.contains(".")} ; then
	rm  TMP/jeter.bam
fi


mv -v TMP FASTQS
"""
}


process BWA_MEM_01 {
tag "${sample} ${R1.name} ${R2.name.contains(".")?R2.name:""}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bwa.yml"
input:
	path(genome)
	path(bwa_index)
	path(bed)
	tuple val(sample),path(R1),path(R2)
output:
	tuple val(sample),path("${sample}.*")
script:
	def fasta = genome.find{it.name.endsWith("a")}
	def ID = sample
	def CN = "Nantes"
	def PL = "ILLUMINA"
	def LB = sample
	def cpus1 = ((task.cpus?:2) as int)
	def cpus2 = cpus1 -1
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP
set -x

BWA_INDEX=`find ${bwa_index}/ -name "*.amb" | sed 's/\\.amb\$//'`

bwa mem -t ${cpus2} \\
	-R '@RG\\tID:${ID}\\tSM:${sample}\\tLB:${LB}\\tCN:${CN}\\tPL:${PL}' \\
	"\${BWA_INDEX}" "${R1}" ${R2.name.contains(".")?"\"${R2}\"":""} |\\
	samtools view -O BAM -o TMP/jeter.bam

if ${R2.name.contains(".")} ; then

	# collate
	samtools collate -l 5 --threads ${task.cpus}  --output-fmt BAM --no-PG --reference "${reference}" -T TMP/tmp.collate -o TMP/jeter2.bam TMP/jeter.bam
	mv TMP/jeter2.bam TMP/jeter.bam

	# fixmate
	samtools fixmate --threads  ${task.cpus} -mc --output-fmt BAM TMP/jeter.bam TMP/jeter2.bam
	mv TMP/jeter2.bam TMP/jeter.bam

fi


samtools sort --threads ${task.cpus} -o TMP/jeter2.bam -O BAM -T TMP/tmp TMP/jeter.bam
mv TMP/jeter2.bam TMP/jeter.bam

samtools index  -@ ${task.cpus} TMP/jeter.bam

mv "TMP/jeter.bam" "${sample}.sorted.bam"
mv "TMP/jeter.bam.bai" "${sample}.sorted.bam.bai"
"""
}

process MERGE_BAMS {
label "process_short"
tag "${sample}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	path(genome)
	tuple val(sample),path("BAMS/*")
output:
	tuple val(sample),path("${sample}.*"),emit:output
script:
	def fasta = genome.find{it.name.endsWith("a")}
"""
mkdir -p TMP
find BAMS -name "*am" > TMP/jeter.list

samtools merge \\
	--threads ${task.cpus} \\
	--reference "${fasta}" \\
	-b TMP/jeter.list \\
	-O  BAM \\
	-o TMP/jeter.bam

# flags with ms
samtools view -u -O BAM  -e '[ms]' TMP/jeter.bam |\\
	samtools markdup \\
	-s --json -f TMP/stats.json \\
	-T TMP/tmp \\
	--threads ${task.cpus} \\
	'-' "TMP/tmp1.bam"

# without ms
samtools view  -O BAM -e '![ms]' -o TMP/tmp2.bam TMP/jeter.bam

samtools merge \\
	--threads ${task.cpus} \\
	--reference "${fasta}" \\
	-O CRAM \\
	--write-index \\
	-o "TMP/${sample}.cram" TMP/tmp1.bam TMP/tmp2.bam

mv TMP/stats.json markdup.${sample}.json
mv -v TMP/${sample}.* ./
"""
}


process COLLATE_AND_MAP_BWA {
label "process_short"
tag "${sample} ${bam.name}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bwa.yml"
input:
	path(genome)
	path(bwa_index)
	path(bed1)
	tuple val(sample),path(bam),path(bai),path(fasta),path(fai),path(bed0)
output:
	tuple val(sample),path("${sample}.sorted.*"),emit:output
script:
	def reference = genome.find{it.name.endsWith("a")}
	def ID = sample
	def CN = "Nantes"
	def PL = "ILLUMINA"
	def LB = sample
	def cpus1 = ((task.cpus?:2) as int)
	def cpus2 = java.lang.Math.max(1,cpus1 -4)
	def fetch_pairs = task.attempt==1?" --fetch-pairs ":""
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP
set -x
BWA_INDEX=`find ${bwa_index}/ -name "*.amb" | sed 's/\\.amb\$//'`


	samtools view -u \\
		 --reference "${fasta}" \\
		-O BAM \\
		${bed0.name.contains(".")?"${fetch_pairs} -M --regions-file \"${bed0}\"":""} \\
		"${bam}" |\\
	samtools collate -f \\
	        -O -u --no-PG \\
        	--reference "${fasta}" \\
	        - TMP/tmp.collate |\\
        samtools fastq -n  -0 TMP/R0.fq.gz -s TMP/Rs.fq.gz - |\\
	bwa mem -t ${cpus2} \\
		-p \\
		-R '@RG\\tID:${ID}\\tSM:${sample}\\tLB:${LB}\\tCN:${CN}\\tPL:${PL}' \\
		"\${BWA_INDEX}" - |\\
	samtools view -O BAM ${bed1.name.contains(".")?"-L \"${bed1}\"":""} -o TMP/jeter.bam


	# collate
	samtools collate -l 5 --threads ${task.cpus}  --output-fmt BAM --no-PG --reference "${reference}" -T TMP/tmp2.collate -o TMP/jeter2.bam TMP/jeter.bam
	mv TMP/jeter2.bam TMP/jeter.bam

	# fixmate
	samtools fixmate --threads  ${task.cpus} -mc --output-fmt BAM TMP/jeter.bam TMP/jeter2.bam
	mv TMP/jeter2.bam TMP/jeter.bam

	# sort
	samtools sort --threads ${task.cpus} -o TMP/jeter2.bam -O BAM -T TMP/tmp TMP/jeter.bam
	mv TMP/jeter2.bam TMP/jeter.bam


	# map single reads
	gunzip -c TMP/R0.fq.gz TMP/Rs.fq.gz |\\
		bwa mem -t ${cpus2} \\
			-R '@RG\\tID:${ID}\\tSM:${sample}\\tLB:${LB}\\tCN:${CN}\\tPL:${PL}' \\
			"\${BWA_INDEX}" - |\\
		samtools view -u -O BAM ${bed1.name.contains(".")?"-L \"${bed1}\"":""} |\\
		samtools fixmate -umc  - - |\\
		samtools sort -o TMP/jeter0.bam -O BAM  -T TMP/tmp3 -


	# merge paired and single end
	samtools merge --threads  ${task.cpus} -o TMP/merged.bam TMP/jeter.bam TMP/jeter0.bam
	mv TMP/merged.bam TMP/jeter.bam


samtools index  -@ ${task.cpus} TMP/jeter.bam

MD5=`cat TMP/jeter.bam | md5sum | cut -d ' ' -f1`

mv "TMP/jeter.bam" "${sample}.sorted.\${MD5}.bam"
mv "TMP/jeter.bam.bai" "${sample}.sorted.\${MD5}.bam.bai"
"""
}

