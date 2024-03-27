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


include {getVersionCmd} from '../../modules/utils/functions.nf'
include {SAMTOOLS_SAMPLES02} from '..//samtools/samtools.samples.02.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {MANTA_SINGLE_01} from  '../../modules/manta/manta.single.01.nf'
include {SIMPLE_ZIP_01} from '../../modules/utils/zip.simple.01.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {TRUVARI_01} from '..//truvari/truvari.01.nf'

workflow MANTA_SINGLE_SV01 {
	take:
		meta
		genomeId
		bams
		bed
		bed_tbi
	main:
		version_ch = Channel.empty()


		all_samples_ch = SAMTOOLS_SAMPLES02([allow_multiple_references:false,with_header:true,allow_duplicate_samples:false],genomeId,bams)
		version_ch= version_ch.mix(all_samples_ch.version)

		each_sample_bam_ch = all_samples_ch.output.splitCsv(header:true,sep:'\t')
	
		
		manta_ch = MANTA_SINGLE_01([:],genomeId,bed,bed_tbi,each_sample_bam_ch)
		version_ch= version_ch.mix(manta_ch.version.first())
	

		all_manta_files_ch = manta_ch.output.
                      map{T->T[1]}.
                      splitCsv(header:true,sep:'\t').
                      map{T->T.vcf}
		
		to_zip = all_manta_files_ch 

		file_list_ch = COLLECT_TO_FILE_01([:],all_manta_files_ch.collect())

		
		if(params.with_merge_manta_vcf==false) {
			merge_vcf = Channel.empty()
			merge_vcf_index = Channel.empty()
			}
		else
			{

			truvari_ch = TRUVARI_01([:] , genomeId, file_list_ch.output)
			to_zip = to_zip.mix(truvari_ch.version)

			merge_vcf = truvari_ch.vcf
			merge_vcf_index = truvari_ch.index

			to_zip = to_zip.mix(merge_vcf).
					mix(merge_vcf_index)
			}


		
		version_merge = MERGE_VERSION("manta",version_ch.collect())
		to_zip = to_zip.mix(version_merge.version)
		
		zip_ch = SIMPLE_ZIP_01(["compression_level":"0"],to_zip.collect())
		
	emit:
		version = version_merge.version
		zip = zip_ch.zip
		merge_vcf = merge_vcf
		merge_vcf_index = merge_vcf_index
		manta_files = file_list_ch.output
	}

