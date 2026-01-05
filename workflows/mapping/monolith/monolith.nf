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


include {runOnComplete;moduleLoad;getVersionCmd} from '../../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.02.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'

/**

map fastq using BWA, all in the same process

*/
workflow {
	ch1_ch  = MONOLITH([:], params.genomeId, Channel.fromPath(params.samplesheet), file(params.bed))
	html = VERSION_TO_HTML(ch1_ch.version)
	}


runOnComplete(workflow);


workflow MONOLITH {
	take:
		meta
		genomeId
		fastqs
		bed
	main:
		version_ch = Channel.empty()

		fastq_ch = fastqs.splitCsv(header:true,sep:'\t').
			map{T->{
				def h = [:]
				if(!T.containsKey("ID") || T.ID.equals(".") || T.ID.trim().isEmpty()) {
					h=h.plus("ID":T.sample)
					}
				if(!T.containsKey("LB") || T.LB.equals(".") || T.LB.trim().isEmpty()) {
					h=h.plus("LB":T.sample)
					}
				if(!T.containsKey("CN") || T.CN.equals(".") || T.CN.trim().isEmpty()) {
					h=h.plus("CN":params.bwa.mem.CN)
					}
				if(!T.containsKey("PL") || T.PL.equals(".") || T.PL.trim().isEmpty()) {
					h=h.plus("PL":params.bwa.mem.PL)
					}
				if((!T.containsKey("R2") || T.R2.equals(".")) && (!T.containsKey("R1") || T.R1.equals(".")) && T.containsKey("ora")) {
					h= h.plus("R1":T.ora);
					h= h.plus("R2":".");
					}
				else if(!T.containsKey("R2") || T.R2.equals(".") || T.R2.trim().isEmpty() || T.R1.endsWith(".ora")) {
					h= h.plus("R2":".");
					}
				return T.plus(h);
				}
			}.map{T->[T.sample,T]}.
			groupTuple();

		bwa_ch = APPLY_BWA(meta, genomeId, bed, fastq_ch)
		version_ch = version_ch.mix(bwa_ch.version)


		version_ch = MERGE_VERSION("Monolith",version_ch.collect())
	emit:
		version = version_ch
	}




process APPLY_BWA {
tag "${sample} ${file(L[0].R1).name} ${L[0].R2.equals(".")?"":file(L[0].R2).name} N=${L.size()}"
afterScript "rm -rf TMP"
// memory set in config
// cpus set in config
input:
	val(meta)
	val(genomeId)
	val(bed)
	tuple val(sample),val(L)
output:
	tuple path("${params.prefix?:""}${sample}.${genomeId}.cram"),path("${params.prefix?:""}${sample}.${genomeId}.cram.crai"),emit:bam
	path("version.xml"),emit:version
script:
	def row = L[0]
	if(!row.R1.endsWith(".ora") && !row.R1.endsWith("q.gz")) throw new IllegalArgumentException("fastq must end with q.gz or .ora")
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
	def bwa_reference = genome.bwa_reference?:reference
	def description = params.description?:""
	def filterbed = file(bed).name.equals("NO_FILE")?"": "${params.htsplusplus.executable} samviewwithmate -B '${bed}' -O BAM --compression 0 |"
"""
hostname 1>&2
${moduleLoad("samtools bwa htsplusplus")}
set -o pipefail
mkdir -p TMP/FASTP
set -x

cat << EOF > TMP/samplesheet.tsv
${L.collect{T->(T.ID?:sample)+"\t"+(T.LB?:sample)+"\t"+ (T.CN?:"BirdNantes")+"\t"+(T.PL?:"ILLUMINA")+"\t"+T.R1+"\t"+T.R2}.join("\n")}
EOF

test -s TMP/samplesheet.tsv

i=1
cat TMP/samplesheet.tsv  | while IFS=\$'\\t' read ID LB CN PL R1 R2
do
	rm -vf TMP/*q.gz 1>&2

	if [[ "\${R1}"  =~ .ora\$  ]] ; then
		${moduleLoad("orad")}
		orad "\${R1}" -V -N -t '${task.cpus}' -P TMP
	else
		ln --verbose -s "\${R1}" TMP/ 1>&2
		if [[ "\${R2}" != "." ]] ; then
			ln --verbose -s "\${R2}" TMP/  1>&2
		fi
	fi

	ls -lah TMP 1>&2
	du -hs TMP 1>&2
	find TMP -name "*q.gz" | grep -F q.gz 1>&2


	if ${params.with_fastp} ; then
		${moduleLoad("fastp")}
		if [[ \$(find TMP  -name "*q.gz"  | wc -l) -eq 1  ]] ; then
			R1=`find TMP  -name "*q.gz"`
			fastp --thread ${task.cpus}  -i "\${R1}" -o TMP/FASTP/fastp.R1.fq.gz
			rm -v "\${R1}" 1>&2
			mv -v TMP/FASTP/fastp.R1.fq.gz "\${R1}" 1>&2
		else
			R1=`find TMP  -name "*q.gz" | sort -V | head -n1`
			R2=`find TMP  -name "*q.gz" | sort -V | tail -n1`
			fastp --thread ${task.cpus} -i "\${R1}" -I "\${R2}" -o TMP/FASTP/fastp.R1.fq.gz -O TMP/FASTP/fastp.R2.fq.gz
			rm -v "\${R1}" 1>&2
			rm -v "\${R2}" 1>&2
			mv -v TMP/FASTP/fastp.R1.fq.gz "\${R1}"	1>&2
			mv -v TMP/FASTP/fastp.R2.fq.gz "\${R2}"	1>&2
		fi
		du -hs TMP 1>&2
	fi

	bwa mem  -t '${task.cpus}' \
		${description.isEmpty()?"":"-H '@CO\\t${description}'"} \
		-R "@RG\\\\tID:\${ID}\\\\tSM:${sample}\\\\tLB:\${LB}\\\\tCN:\${CN}\\\\tPL:\${PL}" \
		"${bwa_reference}" \
		`find TMP  -name "*q.gz"| sort -V  ` | ${filterbed}  \

TODO add collate

		samtools fixmate -m -c -O BAM - TMP/jeter.bam

	rm -vf TMP/*q.gz 1>&2
	du -hs TMP 1>&2

	samtools sort  -m '${task.memory.giga}G' --threads '${task.cpus}' -o TMP/jeter2.bam -O BAM -T TMP/tmp TMP/jeter.bam
	mv -v TMP/jeter2.bam TMP/jeter.bam 1>&2
	du -hs TMP 1>&2

	mv -v TMP/jeter.bam TMP/chunck.\${i}.bam 1>&2
	echo "TMP/chunck.\${i}.bam" >> TMP/bams.list

	i=\$((i + 1))
done

if [[ \$(wc -l < TMP/bams.list) -gt 1  ]] ; then

	samtools merge --threads ${task.cpus} -O BAM -f -o TMP/jeter.bam -b TMP/bams.list
	rm -v TMP/chunck.*.bam 1>&2
	du -hs TMP 1>&2
else
	mv -v  TMP/chunck.1.bam TMP/jeter.bam 1>&2

fi

samtools markdup  --reference "${reference}" --threads ${task.cpus} -T TMP/tmp -O "CRAM,level=9" --write-index  TMP/jeter.bam TMP/jeter2.cram
rm -v TMP/jeter.bam 1>&2
du -hs TMP 1>&2

mv -v "TMP/jeter2.cram" "${params.prefix?:""}${sample}.${genomeId}.cram" 1>&2
mv -v "TMP/jeter2.cram.crai" "${params.prefix?:""}${sample}.${genomeId}.cram.crai" 1>&2
du -hs TMP 1>&2

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">map fastq to reference with bwa</entry>
	<entry key="sample">${sample}</entry>
	<entry key="reference">${reference}</entry>
	<entry key="bwa.version">${getVersionCmd("samtools bwa")}</entry>
</properties>
EOF
"""
}
