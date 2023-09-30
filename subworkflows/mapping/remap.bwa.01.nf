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
include {MERGE_VERSION as MERGE_VERSIONA; MERGE_VERSION as MERGE_VERSIONB; MERGE_VERSION as MERGE_VERSIONC} from '../../modules/version/version.merge.02.nf'
include {moduleLoad;getVersionCmd} from '../../modules/utils/functions.nf'
include {SAMTOOLS_SAMPLES} from '../samtools/samtools.samples.03.nf'
include {SAMTOOLS_COLLATE} from '../../modules/samtools/samtools.collate.01.nf'


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


		remap_ch = REMAP_ONE([:], genomeId, bed, samples_ch.rows.map{T->T.plus("sample":T.new_sample)} )
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

	collate_ch = COLLATE_BAMS([:], bed, row)
	version_ch = version_ch.mix(collate_ch.version)

			
	bam1_ch = MAP_BWA_01([:],genomeId, collate_ch.output)
	version_ch = version_ch.mix(bam1_ch.version)

	version_ch = MERGE_VERSIONB("Remap bam on another reference", version_ch.collect())

emit:
	version=version_ch
	bams = bam1_ch.bams
}


workflow COLLATE_BAMS {
take:
	meta
	bed
	row
main:
	
     	version_ch = Channel.empty()

        collate_ch = SAMTOOLS_COLLATE([:], bed, row)
        version_ch = version_ch.mix(collate_ch.version)


        src1_ch = collate_ch.output.
		map{T->T[0].plus("R1":T[1],"R2":T[2])}.
                filter{T->!(isEmptyGz(T.R1) && isEmptyGz(T.R2))}

        src2_ch = collate_ch.output.
		flatMap{T-> [
			T[0].plus("R1":T[3]),
			T[0].plus("R1":T[4])
			]}.
                filter{T->!isEmptyGz(T.R1)}




        src3_ch =  src1_ch.mix(src2_ch).
                map{T->{
                        T.remove("bam");
                        T.remove("new_sample");
                        T.remove("genomeId");
                        T.remove("fasta");
                        T.remove("gtf");
                        T.remove("reference");
                        return T;
                        }}


	version_ch = MERGE_VERSIONC("Extract FASTQ from BAM", version_ch.collect())
emit:
	output = src3_ch
	version = version_ch

}


