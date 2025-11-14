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

/** convert ORA to FASTQ using illumina ORAD */
workflow ORA_TO_FASTQ {
	take:
		workflow_metadata /* must be provided using a Channel.of(xxx) */
		ora_files /* [meta,oraR0] or [meta,oraR1,oraR2] */
	main:
		
		version_ch = Channel.empty()
		
		meta2 =ora_files
				.count()
				.filter{it>0}
				.map{[:].plus(workflow_metadata)}.first()


		DOWNLOAD_ORAD(meta2)	
		version_ch = version_ch.mix(DOWNLOAD_ORAD.out.versions)

		ch1 = ora_files
			.flatMap{
				if(it.size()==2) return [[it[0].plus("side":"R0"),it[1]]]
				if(it.size()!=3)  {
					log.err("ORA_TO_FASTQ: Expected 2 or 3 values in list: ${it}");
					exit -1
					}
				return [
					[it[0].plus(side:"R1"),it[1]],
					[it[0].plus(sdie:"R2"),it[2]]
					];
				}

		ORAD_INFO(DOWNLOAD_ORAD.out.oradir, ch1)
		version_ch = version_ch.mix(ORAD_INFO.out.versions)
		
		info_ch = ORAD_INFO.out.info
			.splitCsv(header:true,sep:',',strip:true)
			.map{meta,info,ora->[meta.plus(original_file_name:info.Original_file_name),info,ora]}
			.branch{meta,info,ora->
				interleaved: info.Contains_interleaved_data=="YES" && meta.side=="R0"
				single: info.Contains_interleaved_data=="NO" && meta.side=="R0"
				standard : info.Contains_interleaved_data=="NO" && (meta.side=="R1" || meta.side=="R2")
				other: false
			}

		info_ch.other.map{
			log.warn("ORA_TO_FASTQ: Error Contains_interleaved_data = unknown'${it}'.");
			exit -1
			}


		paired_end = Channel.empty()
		single_end = Channel.empty()
		fastqs = Channel.empty()

		/**
		 *
		 * INTERLEAVED DATA
		 *
		 */ 

		RUN_ORAD_INTERLEAVED(
			DOWNLOAD_ORAD.out.oradir ,
			info_ch.interleaved
				.map{meta,info,ora->[[id:meta.id],ora]}
				.groupTuple()
				.map{meta,ora_files->[meta.plus(interleaved:true),ora_files.sort()]}
			)
		version_ch = version_ch.mix(RUN_ORAD_INTERLEAVED.out.versions)


		

		/**
		 *
		 * SINGLE-END DATA
		 *
		 */ 

		RUN_ORAD_SE(
			DOWNLOAD_ORAD.out.oradir ,
			info_ch.interleaved
				.map{meta,info,ora->[[id:meta.id],ora]}
				.groupTuple()
				.map{meta,ora_files->[meta.plus(interleaved:false),ora_files.sort()]}
			)
		version_ch = version_ch.mix(RUN_ORAD_SE.out.versions)

		/**
		 *
		 * PAIRED-END DATA DATA
		 *
		 */ 
		RUN_ORAD_PE(
			info_ch.interleaved
				.map{meta,info,ora->[[id:meta.id,side:meta.side],ora]}
				.groupTuple()
				.map{meta,ora_files->[meta.plus(interleaved:false),ora_files.sort()]}
			)
		version_ch = version_ch.mix(RUN_ORAD_PE.out.versions)


		paired_R1 = RUN_ORAD_PE.out.fastq
				.filter{meta,_fastq->meta.side=="R1"}
				.map{meta,fastq->[meta.findAll{k,v->!k.matches("(side|interleaved)")},fastq]}
		
		paired_R2 = v.out.fastq
				.filter{meta,_fastq->meta.side=="R2"}
				.map{meta,fastq->[meta.findAll{k,v->!k.matches("(side|interleaved)")},fastq]}


		paired_ch = paired_R1.join(paired_R2,failOnMismatch:true,failOnDuplicate:true)


		if(1==2) {
			RUN_ORAD( DOWNLOAD_ORAD.out.oradir ,ora_files)
			version_ch = version_ch.mix(RUN_ORAD.out.versions)

		
			ch1 = RUN_ORAD.out.fastq.branch{v->
				paired_end: v[1].size()==2
				single_end: v[1].size()==1
				other: true
				}		

			ch1.other.map{
				throw new IllegalStateException("orad unexpected output : ${it}.");
				}

			paired_end = ch1.paired_end.map{
				def L = it[1].sort() // sort R1 and R2
				return [it[0],L[0],L[1]];
				}
			single_end = ch1.single_end.map {
				return [it[0],it[1][0]];
				}
			fastqs = paired_end.mix( single_end.map{m,fq->[m,fq,[] ]} )
			}
	emit:
		versions = version_ch
		paired_end
		single_end
		fastqs
	}


