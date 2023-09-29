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


include {MAP_BWA_01} from '../../subworkflows/mapping/map.bwa.01.nf'
include {MERGE_VERSION as MERGE_VERSIONA; MERGE_VERSION as MERGE_VERSIONB} from '../../modules/version/version.merge.02.nf'
include {isBlank;moduleLoad;getVersionCmd} from '../../modules/utils/functions.nf'
include {SAMTOOLS_SAMPLES} from '../samtools/samtools.samples.03.nf'


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
	}

workflow REMAP_BWA_01 {
	take:
		meta
		genomeId
		bams
		bed
	main:
		version_ch = Channel.empty()


		samples_ch = SAMTOOLS_SAMPLES([:],bams)
		version_ch = version_ch.mix(samples_ch.version)


		remap_ch = REMAP_ONE([:], genomeId, bed, samples_ch.rows )
		version_ch = version_ch.mix(remap_ch.version)
		version_ch = MERGE_VERSIONA("Extract Fastq from bam and remap", version_ch.collect())

	emit:
		version = version_ch
		bams = remap_ch.bams 
	}

workflow REMAP_ONE {
take:
	meta
	genomeId
	bed
	row
main:
	version_ch = Channel.empty()

	collate_ch = COLLATE_AND_FASTQ([:], bed, row)
	version_ch = version_ch.mix(collate_ch.version)


	src1_ch = collate_ch.paired.splitCsv(header:true,sep:'\t').
		filter{T->!(isEmptyGz(T.R1) && isEmptyGz(T.R2))}

	src2_ch = collate_ch.single.splitCsv(header:true,sep:'\t').
		filter{T->!isEmptyGz(T.R1)}

	src3_ch =  src1_ch.mix(src2_ch).
		map{T->row.plus(T)}.		
		map{T->{
                	T.remove("bam");
                	T.remove("genomeId");
	                T.remove("fasta");
	                T.remove("gtf");
	                T.remove("reference");
        	        return T;
                	}}.
	        map{T->T.plus("genomeId":genomeId)}

/*
			
	bam1_ch = MAP_BWA_01([:],genomeId, src3_ch)

	version_ch = version_ch.mix(bam1_ch.version)

	version_ch = MERGE_VERSIONB("Remap bam on another reference", version_ch.collect())
emit:
	version= version_ch.version
	bams = bam1_ch.bams
*/
	bamsxx  = Channel.empty()
emit:
	version=version_ch
	bams = bamsxx
}



process COLLATE_AND_FASTQ {
tag "${row.new_sample} ${file(row.bam).name}"
afterScript "rm -rf TMP TMP2"
memory "5g"
cpus 4
input:
	val(meta)
	path(bed)
	val(row)
output:
	path("FASTQS/${row.new_sample}.paired.fastq.tsv"),emit:paired
	path("FASTQS/${row.new_sample}.single.fastq.tsv"),emit:single
	path("version.xml"),emit:version
script:
	def sample = row.new_sample
"""
hostname 1>&2
set -o pipefail
${moduleLoad("samtools/")}

mkdir -p TMP TMP2


# extract reads
if ${!bed.name.equals("NO_FILE")} ; then
	samtools view  --threads ${task.cpus} --reference "${row.reference}" -M -O BAM -o TMP/jeter.bam --regions-file "${bed}" "${row.bam}"
fi

# show size
ls -lah  "${!bed.name.equals("NO_FILE") ? "TMP/jeter.bam":"${row.bam}"}" 1>&2


# collate
samtools collate -f --threads ${(task.cpus as int) -1} -O -u --no-PG --reference "${row.reference}" "${!bed.name.equals("NO_FILE") ? "TMP/jeter.bam":"${row.bam}"}" TMP/tmp.collate |\
	samtools fastq -n --threads 1 -1 TMP/jeter.paired.R1.fq.gz -2 TMP/jeter.paired.R2.fq.gz -s "TMP/jeter.singleton.fq.gz" -0 "TMP/jeter.other.fq.gz"


# merge both singleton or other

cat TMP/jeter.singleton.fq.gz TMP/jeter.other.fq.gz > TMP/jeter.gz
mv -v TMP/jeter.gz TMP/jeter.singleton.fq.gz

mv -v TMP2 FASTQS

find \${PWD}/FASTQS -type f -name "${sample}.paired*.R1.fastq.gz" -o -name "${sample}.paired*.R2.fastq.gz" | LC_ALL=C sort -T TMP | paste - - |\
	awk -F '\t' 'BEGIN{printf("R1\tR2\\n");} {printf("%s\t%s\\n",\$1,\$2);}' > FASTQS/${sample}.paired.fastq.tsv

find \${PWD}/FASTQS -type f -name "${sample}.singleton*.fastq.gz" | \
	awk -F '\t' 'BEGIN{printf("R1\\n");} {printf("%s\\n",\$0);}' > FASTQS/${sample}.single.fastq.tsv

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

