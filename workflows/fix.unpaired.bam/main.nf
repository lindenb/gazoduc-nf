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


include {dumpParams;getVersionCmd;runOnComplete;moduleLoad;isBlank} from '../../modules/utils/functions.nf'


if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}

String bamPrefix(String fn) {
	fn = file(fn).name;
	if(fn.contains("/")) throw new IllegalArgumentException(fn);
	if(fn.endsWith(".bam")) fn = fn.substring(0,fn.length()-4);
	if(fn.endsWith(".cram")) fn = fn.substring(0,fn.length()-5);
	if(fn.endsWith(".hs38me")) fn =  fn.substring(0,fn.length()-7);
	if(!fn.endsWith(".")) fn+=".";
	return fn;
	}

workflow {
	ch = FIX_UNPAIRED_BAMS([:],params.fasta, file(params.bams) )
	}

runOnComplete(workflow);

workflow FIX_UNPAIRED_BAMS {
    take:
	    meta
	    fasta
	    bams
    main:

		ch0 = Channel.fromPath(bams).splitText().map{it.trim()};
		ch1  = BAM_WITH_SINGLE(fasta, ch0 )
		ch4 = FLAGSTATS(fasta,ch0)

		ch2 = ch1.output.
			splitCsv(header:false, sep: '\t').
			map{T->[
				sample: T[0],
				bam: T[1],
				fasta : fasta
				]}


		ch3= FIX_BAM(ch2)

    }


/** scan bam with single end */
process BAM_WITH_SINGLE {
tag "${file(bam).name}"
label  "process_medium"
cpus 6
input:
	val(fasta)
	val(bam)
output:
	path("samples.tsv"),emit:output
script:
"""
module load samtools
samtools view  --threads ${task.cpus}  -F 1 -T "${fasta}" "${bam}" |\\
	awk '{printf("${bam}\\n"); exit(0)}' |\\
	samtools samples > samples.tsv
"""
}


process FLAGSTATS {
tag "${file(bam).name}"
label  "process_medium"
cpus 6
input:
        val(fasta)
        val(bam)
output:
        path("${bamPrefix(bam)}flags.tsv"),emit:output
script:
"""
module load samtools
samtools flagstats -O tsv --threads ${task.cpus} "${bam}" |\\
        awk -F '\t' '{printf("%s\t${bam}\\n",\$0);}' > "${bamPrefix(bam)}flags.tsv"
"""
}

process FIX_BAM {
tag "${row.sample} ${file(row.bam).name}"
afterScript "rm -rf TMP"
label  "process_medium"
cpus 10
input:
	val(row)
output:
	tuple val(row),path("${bamPrefix(row.bam)}cram"),path("${bamPrefix(row.bam)}cram.crai"),emit:output
when:
	false
script:
	def overflow_size = 600_000
"""
hostname 1>&2
module load samtools bwa sambamba
mkdir -p TMP
set -x

# extract single and paired end
samtools view --remove-flags 1024 --reference "${fasta}" --threads ${task.cpus} --require-flags 1 -O BAM --output TMP/paired.bam  ----output-unselected TMP/single.bam "${row.bam}"

# sort single
samtools collate --threads ${task.cpus} --output-fmt BAM -o TMP/jeter.bam -f -O -u --no-PG --reference "${fasta}" TMP/single.bam TMP/tmp.collate 
mv TMP/jeter.bam TMP/single.bam

# map single to fastq
samtools fastq -N -1 TMP/jeter.R1.fq.gz -2 TMP/jeter.R2.fq.gz -s TMP/jeter.Rx.fq.gz -0 TMP/jeter.R0.fq.gz -n TMP/single.bam

cat TMP/jeter.Rx.fq.gz TMP/jeter.R0.fq.gz | gunzip -c | paste - - - - | LC_ALL=C sort -t '\t' -T TMP -k1,1 > TMP/jeter.fq.tsv
wc -l TMP/jeter.fq.tsv 1>&2
head -n 40 TMP/jeter.fq.tsv 1>&2


# join, create a interleaved fastq
LC_ALL=C join -t '\t' -1 1 -2 1 -o '1.1,1.2,1.3,1.4,2.1,2.2,2.3,2.4' TMP/jeter.fq.tsv TMP/jeter.fq.tsv |\\
	awk '(\$1==\$5 && \$2 < \$6)' |\
	tr "\t" "\\n" > TMP/jeter.fastq

wc -l TMP/jeter.fastq 1>&2
head -n 40 TMP/jeter.fastq 1>&2

bwa mem -t '${task.cpus}' -p -R '@RG\\tID:${row.sample}\\tSM:${row.sample}\\tLB:${row.sample}\\tCN:NantesBird\\tPL:ILLUMINA' "${fasta}"  TMP/jeter.fastq |\\
	samtools view -O BAM -o TMP/jeter.bam
	
samtools sort -T TMP/sort  -@ ${task.cpus} -O BAM --reference "${fasta}" -o jeter2.bam   TMP/jeter.bam
mv -v jeter2.bam TMP/jeter.bam
mv -v TMP/jeter.bam TMP/single.bam


samtools merge --threads ${task.cpus}  --reference "${fasta}" -O BAM -o TMP/jeter.bam  TMP/paired.bam TMP/single.bam
rm -vf TMP/paired.bam TMP/single.bam


# avoid too many files open. ulimit doesn't work...
ulimit -s unlimited 

# markdup
sambamba markdup --overflow-list-size ${overflow_size} --tmpdir=TMP -t ${task.cpus} TMP/jeter.bam TMP/jeter2.bam
mv -v TMP/jeter2.bam TMP/jeter.bam


samtools view --write-index --reference "${fasta}" -O "CRAM,level=9" -o TMP/jeter.cram TMP/jeter.bam

samtools flagstat --threads ${task.cpus} TMP/jeter.cram > "${bamPrefix(row.bam)}flagstats.txt"

mv TMP/jeter.cram "${bamPrefix(row.bam)}cram"
mv TMP/jeter.cram.crai "${bamPrefix(row.bam)}cram.crai"
"""
}
