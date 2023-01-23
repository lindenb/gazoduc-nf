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


def gazoduc = gazoduc.Gazoduc.getInstance()

gazoduc.build("split_fastqs_count", 100).
	desc("Split a FASTQ produced by samtools fastq into 'N' files").
	setInt().
	menu("fastq").
	put()

gazoduc.build("split_fastq_ignore_if_size", 10000000L).
	desc("Do NOT split the fastq file if it's compressed size if lower than 'x' bytes.").
	setLong().
	menu("fastq").
	put()

include {MAP_BWA_01} from '../../subworkflows/mapping/map.bwa.01.nf'
include {COMPILE_FASTQ_SPLIT2FILE} from '../../modules/fastq/fastq.splitter.02.nf'
include {MERGE_VERSION as MERGE_VERSIONA; MERGE_VERSION as MERGE_VERSIONB} from '../../modules/version/version.merge.nf'
include {isBlank;moduleLoad;getVersionCmd} from '../../modules/utils/functions.nf'
include {SAMTOOLS_SAMPLES_01} from '../samtools/samtools.samples.01.nf'

/*
boolean isEmptyGz(Object filename) {
	final java.nio.file.Path p;
	if(filename instanceof java.nio.file.Path) {
		p = (java.nio.file.Path)filename;
		}
	else
		{
		p = java.nio.file.Paths.get(filename.toString());
		}
	try (java.util.zip.GZIPInputStream in = new  java.util.zip.GZIPInputStream(java.nio.file.Files.newInputStream(p)) ) {
		return in.read()==-1;
		}
	}*/

workflow REMAP_BWA_01 {
	take:
		meta
		reference_in
		reference_out
		bams
		bed
	main:
		version_ch = Channel.empty()

		splitter_ch = COMPILE_FASTQ_SPLIT2FILE(meta)
		version_ch = version_ch.mix(splitter_ch.version)

		samples_ch = SAMTOOLS_SAMPLES_01(
				meta.plus("with_header":true),
				reference_in,
				file("NO_FILE"),
				bams
				)
		version_ch = version_ch.mix(samples_ch.version)

		each_sn_bam =  samples_ch.output.splitCsv(header:true,sep:'\t')

		remap_ch = REMAP_ONE(meta, splitter_ch.executable, reference_out, bed, each_sn_bam )
		version_ch = version_ch.mix(remap_ch.version)
		version_ch = MERGE_VERSIONA(meta, "Remap", "Extract Fastq from bam and remap", version_ch.collect())

	emit:
		version = version_ch
		bams = remap_ch.bams
	}

workflow REMAP_ONE {
take:
	meta
	splitter_executable
	reference_out
	bed
	row
main:
	version_ch = Channel.empty()

	collate_ch = COLLATE_AND_FASTQ(meta, splitter_executable, bed, row)
	version_ch = version_ch.mix(collate_ch.version)


	src1_ch = collate_ch.output.splitCsv(header:true,sep:'\t').map{T->[
		"sample": T.sample,
		"R1": T.R1,
		"R2": T.R2,
		"reference" : reference_out
		]}

	src2_ch = collate_ch.single.splitCsv(header:true,sep:'\t').map{T->[
		"sample": T.sample,
		"R1": T.fastq,
		"reference" : reference_out
		]}
		
	src3_ch = collate_ch.unpaired.splitCsv(header:true,sep:'\t').map{T->[
		"sample": T.sample,
		"R1": T.fastq,
		"reference" : reference_out
		]}
	
	bam1_ch = MAP_BWA_01(meta,reference_out,src1_ch.
		mix(src2_ch).
		mix(src3_ch))
	version_ch = version_ch.mix(bam1_ch.version)

	version_ch = MERGE_VERSIONB(meta, "Remap", "Remap bam on another reference", version_ch.collect())
emit:
	version= version_ch.version
	bams = bam1_ch.bams
}



process COLLATE_AND_FASTQ {
tag "${row.new_sample} ${file(row.bam).name}"
afterScript "rm -rf TMP TMP2"
memory "5g"
cpus 5
input:
	val(meta)
	path(splitter)
	path(bed)
	val(row)
output:
	path("FASTQS/${row.new_sample}.paired.fastq.tsv"),emit:output
	path("FASTQS/${row.new_sample}.singleton.fastq.tsv"),emit:single
	path("FASTQS/${row.new_sample}.other.fastq.tsv"),emit:unpaired
	path("version.xml"),emit:version
script:
	def sample = row.new_sample
"""
hostname 1>&2
set -o pipefail
${moduleLoad("samtools/1.15.1")}

mkdir -p TMP TMP2

# extract reads
if ${!bed.name.equals("NO_FILE")} ; then
	samtools view  --threads ${task.cpus} --reference "${row.reference}" -M -O BAM -o TMP/jeter.bam --regions-file "${bed}" "${row.bam}"
fi

# collate
samtools collate -f --threads 4 -O -u --no-PG --reference "${row.reference}" "${!bed.name.equals("NO_FILE") ? "TMP/jeter.bam":"${row.bam}"}" TMP/tmp.collate |\
	samtools fastq -n --threads 1 -1 TMP/jeter.paired.R1.fq.gz -2 TMP/jeter.paired.R2.fq.gz -s "TMP/jeter.singleton.fq.gz" -0 "TMP/jeter.other.fq.gz"

# split

if [ `stat -c '%s' "TMP/jeter.paired.R1.fq.gz"` -le "${meta.split_fastq_ignore_if_size}" ]
then
	mv -v TMP/jeter.paired.R1.fq.gz "TMP2/${sample}.paired.001.R1.fastq.gz"
	mv -v TMP/jeter.paired.R2.fq.gz "TMP2/${sample}.paired.001.R2.fastq.gz"
else
	${splitter.toRealPath()} -n ${meta.split_fastqs_count} -o "TMP2/${sample}.paired" TMP/jeter.paired.R1.fq.gz  TMP/jeter.paired.R2.fq.gz
fi


if [ `stat -c '%s' "TMP/jeter.singleton.fq.gz"` -le "${meta.split_fastq_ignore_if_size}" ]
then
	mv -v TMP/jeter.singleton.fq.gz "TMP2/${sample}.singleton.001.fastq.gz"
else
	${splitter.toRealPath()} -s -n ${meta.split_fastqs_count} -o "TMP2/${sample}.singleton" TMP/jeter.singleton.fq.gz
fi

if [ `stat -c '%s' "TMP/jeter.other.fq.gz"` -le "${meta.split_fastq_ignore_if_size}" ]
then
	mv -v TMP/jeter.other.fq.gz "TMP2/${sample}.other.001.fastq.gz"
else
	${splitter.toRealPath()} -s -n ${meta.split_fastqs_count} -o "TMP2/${sample}.other" TMP/jeter.other.fq.gz
fi

mv -v TMP2 FASTQS


find \${PWD}/FASTQS -type f -name "${sample}.paired*.R1.fastq.gz" -o -name "${sample}.paired*.R2.fastq.gz" | LC_ALL=C sort -T TMP | paste - - |\
	awk -F '\t' 'BEGIN{printf("sample\tbam\treference\tR1\tR2\\n");} {printf("${sample}\t${row.bam}\t${row.reference}\t%s\t%s\\n",\$1,\$2);}' > FASTQS/${sample}.paired.fastq.tsv

find \${PWD}/FASTQS -type f -name "${sample}.singleton*.fastq.gz" | \
	awk -F '\t' 'BEGIN{printf("sample\tbam\treference\tfastq\\n");} {printf("${sample}\t${row.bam}\t${row.reference}\t%s\\n",\$0);}' > FASTQS/${sample}.singleton.fastq.tsv

find \${PWD}/FASTQS -type f -name "${sample}.other*.fastq.gz" | \
	awk -F '\t' 'BEGIN{printf("sample\tbam\treference\tfastq\\n");} {printf("${sample}\t${row.bam}\t${row.reference}\t%s\\n",\$0);}' > FASTQS/${sample}.other.fastq.tsv



##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">convert BAM to fastq</entry>
	<entry key="sample">${sample}</entry>
	<entry key="bam">${row.bam}</entry>
        <entry key="versions">${getVersionCmd("samtools awk")}</entry>
	<entry key="bed">${bed}</entry>
</properties>
EOF
"""
}

