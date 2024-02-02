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
nextflow.enable.dsl=2


include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {dumpParams;getVersionCmd;runOnComplete;moduleLoad;isBlank} from '../../modules/utils/functions.nf'
include {SAMTOOLS_SAMPLES} from '../../subworkflows/samtools/samtools.samples.03.nf'


if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}


workflow {
	ch = FIX_UNPAIRED_BAMS([:], file(params.bams) )
	html = VERSION_TO_HTML(ch.version)
	}

runOnComplete(workflow);

workflow FIX_UNPAIRED_BAMS {
    take:
	    meta
	    bams
    main:
		version_ch = Channel.empty()

		ch1 = SAMTOOLS_SAMPLES([:], vcf)
		version_ch = version_ch.mix(ch1.version)


		ch4= FIX_BAM([:], ch1.rows)

    emit:
	    version = version_ch
    }


process FIX_BAM {
tag "${row.sample}"
afterScript "rm -rf TMP"
cpus 10
input:
	val(meta)
	val(row)
output:
	path("join.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def genome = params.genomes[row.genomeId]
	def fasta = genome.fasta
"""
hostname 1>&2
${moduleLoad("samtools bwa")}
mkdir -p TMP
set -x

# extract single and paired end
samtools view --reference "${fasta}" --threads ${task.cpus} --require-flags 1 -O BAM --output TMP/paired.bam  ----output-unselected TMP/single.bam "${row.bam}"

# sort single
samtools collate --threads ${task.cpus} --output-fmt BAM -o TMP/jeter.bam -f -O -u --no-PG --reference "${fasta}" TMP/single.bam TMP/tmp.collate 
mv TMP/jeter.bam TMP/single.bam

# map single to fastq
samtools fastq -N -1 TMP/jeter.R1.fq.gz -2 TMP/jeter.R2.fq.gz -s TMP/jeter.Rx.fq.gz -0 TMP/jeter.R0.fq.gz -n TMP/single.bam

cat TMP/jeter.Rx.fq.gz TMP/jeter.R0.fq.gz | gunzip -c | paste - - - - | LC_ALL=C sort -t '\t' -T TMP -k1,1 > TMP/jeter.fq.tsv
wc -l TMP/jeter.fq.tsv 1>&2
head -n 40 TMP/jeter.fq.tsv 1>&2


# join, create a interleaved fastq
LC_ALL=C join -t '\t' -1 1 -2 1 -o '1.1,1.2,1.3,1.4,2.1,2.2,2.3,2.4' TMP/jeter.fq.tsv TMP/jeter.fq.tsv |\
	awk '(\$1==\$5 && \$2 < \$6)' |\
	tr "\t" "\\n" > TMP/jeter.fastq

wc -l TMP/jeter.fastq 1>&2
head -n 40 TMP/jeter.fastq 1>&2

bwa mem -t '${task.cpus}' -p -R '@RG\\tID:${row.sample}\\tSM:${row.sample}\\tLB:${row.sample}\\tCN:NantesBird\\tPL:ILLUMINA' "${fasta}"  TMP/jeter.fastq |\
	samtools view -O BAM -o TMP/jeter.bam
	
samtools sort -T TMP/sort  -@ ${task.cpus} -O BAM --reference "${fasta}" -o jeter2.bam   TMP/jeter.bam
mv jeter2.bam TMP/jeter.bam

samtools fixmate -m -c -O BAM  TMP/jeter.bam jeter2.bam 
mv jeter2.bam TMP/jeter.bam

samtools markdup  --reference "${fasta}"  -T TMP/tmp -O "BAM" --write-index  TMP/jeter.bam TMP/jeter2.bam
mv jeter2.bam TMP/jeter.bam
mv TMP/jeter.bam TMP/single.bam

samtools merge --threads ${task.cpus} --write-index --reference ${fasta} -O "CRAM,level=9" -o TMP/jeter.cram  TMP/paired.bam TMP/single.bam |\
rm -f TMP/paired.bam TMP/single.bam


samtools flagstat --threads ${task.cpus} TMP/jeter.cram > flagsstats.txt

mv TMP/jeter.cram "./${params.prefix?:""}${row.genomeId}.${row.sample}.cram"
mv TMP/jeter.cram.crai "./${params.prefix?:""}${row.genomeId}.${row.sample}.cram.crai"

cat << EOF | paste - - | awk '{printf("mv -v \\"%s\\" \\"%s\\"\\n",\$1,\$2);}' > output.sh
./${params.prefix?:""}${row.genomeId}.${row.sample}.cram
${row.bam}
./${params.prefix?:""}${row.genomeId}.${row.sample}.cram.crai
${row.bam}.crai
EOF



"""
}
