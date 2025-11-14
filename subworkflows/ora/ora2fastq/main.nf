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
include {DOWNLOAD_ORAD                        } from '../../../modules/orad/download'
include {RUN_ORAD  as   RUN_ORAD_INTERLEAVED  } from '../../../modules/orad/run'
include {RUN_ORAD  as   RUN_ORAD_PE           } from '../../../modules/orad/run'
include {RUN_ORAD  as   RUN_ORAD_SE           } from '../../../modules/orad/run'
include {ORAD_INFO                            } from '../../../modules/orad/info'
include {parseBoolean                         } from '../../../modules/utils/functions'
include {isBlank                              } from '../../../modules/utils/functions'
include {verify                               } from '../../../modules/utils/functions'

/** convert ORA to FASTQ using illumina ORAD */
workflow ORA_TO_FASTQ {
	take:
		workflow_metadata /* must be provided using a Channel.of(xxx) */
		ora_files /* [meta,oraR0] or [meta,oraR1,oraR2] */
	main:
		if(workflow_metadata.merge_fastqs) {
			log.warn("ORA_TO_FASTQ: merge_fastqs undefined")
			}
		
		version_ch = Channel.empty()
		
		paired_end = Channel.empty()
		single_end = Channel.empty()
		fastqs = Channel.empty()

		
		meta2 =ora_files
				.count()
				.filter{it>0}
				.map{[:].plus(workflow_metadata)}.first()
	

		/**
		 *
		 * DOWNLOAD ORAD 
		 *
		 */
		DOWNLOAD_ORAD(meta2)	
		version_ch = version_ch.mix(DOWNLOAD_ORAD.out.versions)

		/**
		 *
		 * Add info about side, is it R1 or R2 (or R0 for other)
		 *
		 */
		ch1 = ora_files
			.flatMap{
				verify( isBlank(it[0].side), "side already defined in ${it}.");
				verify( isBlank(it[0].interleaved) , "interleaved already defined in ${it}.");
				if(it.size()==2) return [[it[0].plus("side":"R0"),it[1]]]
				verify( it.size()==3 , "ORA_TO_FASTQ: Expected 2 or 3 values in list: ${it}");
				return [
					[it[0].plus(side:"R1"),it[1]],
					[it[0].plus(side:"R2"),it[2]]
					];
				}
/**
		 *
		 * EXTRACT INFO from ORA file
		 *
		 */
		ORAD_INFO(DOWNLOAD_ORAD.out.oradir, ch1)
		version_ch = version_ch.mix(ORAD_INFO.out.versions)
		
		info_ch = ORAD_INFO.out.info
			.splitCsv(header:true,sep:',',strip:true)
			.map{meta,info,ora->[meta,info,ora]}
			.branch{meta,info,ora->
				boum: isBlank(info.Contains_interleaved_data) || isBlank(meta.side)
				interleaved: info.Contains_interleaved_data=="YES" && meta.side=="R0"
				single: info.Contains_interleaved_data=="NO" && meta.side=="R0"
				paired : info.Contains_interleaved_data=="NO" && (meta.side=="R1" || meta.side=="R2")
				other: false
				}



		
		info_ch.other.mix(info_ch.boum).map{
			verify(false,"ORA_TO_FASTQ: Error Contains_interleaved_data = unknown'${it}'.");
			return -1;
			}


		/**
		 *
		 * INTERLEAVED DATA
		 *
		 */
		if(workflow_metadata.merge_fastqs==null || parseBoolean(workflow_metadata.merge_fastqs)) {
			ch1 = info_ch.interleaved
				.map{meta,_info,ora->[meta.id,meta,ora]}
				.groupTuple()
				.map{_meta_id,metas,files->[metas.sort()[0], files.sort()]}
			}
		else
			{
			ch1 = info_ch.interleaved
				.map{meta,_info,ora->[meta,[ora]]}
			}

		RUN_ORAD_INTERLEAVED(
			DOWNLOAD_ORAD.out.oradir ,
			ch1
			)
		version_ch = version_ch.mix(RUN_ORAD_INTERLEAVED.out.versions)
		paired_end = paired_end.mix(RUN_ORAD_INTERLEAVED.out.fastq
			.map{meta,fastqs->{
				verify(fastqs.size()==2 , "expected two fastq files");
				def L = fastqs.sort();
				return [meta,L[0],L[1]];
				}})
		// paired_end.view()

		/**
		 *
		 * SINGLE-END DATA
		 *
		 */ 
		if(workflow_metadata.merge_fastqs==null || parseBoolean(workflow_metadata.merge_fastqs)) {
			ch2 = info_ch.single
				.map{meta,_info,ora->[meta.id,meta,ora]}
				.groupTuple()
				.map{_meta_id,metas,files->[metas.sort()[0], files.sort()]}
			}
		else
			{
			ch2 = info_ch.single
				.map{meta,_info,ora->[meta,[ora]]}
			}
		RUN_ORAD_SE(
			DOWNLOAD_ORAD.out.oradir ,
			ch2
			)
		version_ch = version_ch.mix(RUN_ORAD_SE.out.versions)
		single_end = single_end.mix(RUN_ORAD_SE.out.fastq
			.map{meta,fastqs->{
				verify( fastqs instanceof List, "should be a List");
				verify( fastqs.size()==1 ,"expected one fastq file");
				return [meta,L[0]];
				}})
		


		/**
		 *
		 * PAIRED-END DATA DATA
		 *
		 */ 
		if(workflow_metadata.merge_fastqs==null || parseBoolean(workflow_metadata.merge_fastqs)) {
			ch3 = info_ch.paired
				.map{meta,_info,ora->[[id:meta.id,side:meta.side],meta,ora]}
				.groupTuple()
				.map{_meta_id_side,metas,files->[meta.sort()[0],files.sort()]}
			}
		else
			{
			ch3 = info_ch.paired
				.map{meta,_info,ora->[meta,[ora]]}
			}

		RUN_ORAD_PE(
			DOWNLOAD_ORAD.out.oradir ,
			ch3
			)
		version_ch = version_ch.mix(RUN_ORAD_PE.out.versions)


		paired_R1 = RUN_ORAD_PE.out.fastq
				.filter{meta,_fastq->meta.side=="R1"}
				.map{meta,fastq->[meta.findAll{k,v->!k.matches("(side)")},fastq]}
		
		paired_R2 = RUN_ORAD_PE.out.fastq
				.filter{meta,_fastq->meta.side=="R2"}
				.map{meta,fastq->[meta.findAll{k,v->!k.matches("(side)")},fastq]}

		paired_ch = paired_R1.join(paired_R2,failOnMismatch:true,failOnDuplicate:true)


		paired_end = paired_end.mix(paired_ch)
				.map{meta,R1,R2->[meta.findAll{k,v->!k.matches("(side)")},R1,R2]}

		single_end = single_end
				.map{meta,R0->[meta.findAll{k,v->!k.matches("(side)")},R0]}

		fastqs = paired_end.mix( single_end.map{m,fq->[m,fq,[] ]} )

	emit:
		versions = version_ch
		paired_end
		single_end
		fastqs
	}


