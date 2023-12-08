/*

Copyright (c) 2023 Pierre Lindenbaum

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
include {moduleLoad;getVersionCmd} from '../utils/functions.nf'



process SAMTOOLS_COLLATE {
tag "${row.sample} ${file(row.bam).name}"
afterScript "rm -rf TMP"
errorStrategy "retry"
maxRetries 2
input:
	val(meta)
	path(bed)
	val(row)
output:
	tuple val(row),
			path("FASTQS/${row.sample}.paired.R1.fq.gz"),
			path("FASTQS/${row.sample}.paired.R2.fq.gz"),
			path("FASTQS/${row.sample}.singleton.fq.gz"),
			path("FASTQS/${row.sample}.other.fq.gz"),emit:output
	path("version.xml"),emit:version
script:
	if(!(row.containsKey("reference") || row.containeKey("genomeId"))) throw new IllegalArgumentException("missing ref/genomeId");
	if(!row.containsKey("sample")) throw new IllegalArgumentException("sample");

	// bug in samtools ?
	def fetch_pairs = task.attempt==1?" --fetch-pairs ":""
	def reference = row.reference?:params.genomes[row.genomeId].fasta
	def level = task.ext.compression_level?:1
	def sample = row.sample
"""
hostname 1>&2
set -o pipefail
${moduleLoad("samtools")}

mkdir -p TMP


# extract reads
if ${!bed.name.equals("NO_FILE")} ; then
	samtools view  --fast --threads ${task.cpus} ${fetch_pairs} --reference "${reference}" -M -O BAM -o TMP/jeter.bam --regions-file "${bed}" "${row.bam}"
fi

# show size
ls -lah  "${!bed.name.equals("NO_FILE") ? "TMP/jeter.bam":"${row.bam}"}" 1>&2


# collate
samtools collate -l ${level} -f --threads ${(task.cpus as int) -1} -O -u --no-PG --reference "${reference}" "${!bed.name.equals("NO_FILE") ? "TMP/jeter.bam":"${row.bam}"}" TMP/tmp.collate |\
	samtools fastq -n --threads 1 -1 TMP/${sample}.paired.R1.fq.gz -2 TMP/${sample}.paired.R2.fq.gz -s "TMP/${sample}.singleton.fq.gz" -0 "TMP/${sample}.other.fq.gz"


if ${!bed.name.equals("NO_FILE")} ; then
	rm  TMP/jeter.bam
fi


mv -v TMP FASTQS

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">convert BAM to fastq</entry>
	<entry key="sample">${sample}</entry>
	<entry key="bam">${row.bam}</entry>
        <entry key="versions">${getVersionCmd("samtools")}</entry>
	<entry key="bed">${bed}</entry>
</properties>
EOF
"""
}

