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
include {DOWNLOAD_ORAD  } from '../../../modules/orad/download'
include {RUN_ORAD       } from '../../../modules/orad/run'

/** convert ORA to FASTQ using illumina ORAD */
workflow ORA_TO_FASTQ {
	take:
		meta /* must be provided using a Channel.of(xxx) */
		ora_files /* [meta,ora] */
	main:
		
		version_ch = Channel.empty()
		
		meta2 =ora_files.count()
				.filter{it>0}
				.combine(meta)
				.map{it[1]}


		DOWNLOAD_ORAD(meta2)	
		version_ch = version_ch.mix(DOWNLOAD_ORAD.out.versions)


		RUN_ORAD( DOWNLOAD_ORAD.out.oradir ,ora_files)
		version_ch = version_ch.mix(RUN_ORAD.out.versions)

	
		ch1 = RUN_ORAD.out.fastq.branch{v->
			paired_end: v[1].size()==2
			single_end: v[1].size()==1
			other: true
			}		

		ch1.other.map{
			throw new IllegalStateException("orad unexpected output : "+it);
			}

		paired_end = ch1.paired_end.map{
			def L = it[1].sort() // sort R1 and R2
			return [it[0],L[0],L[1]];
			}
		single_end = ch1.single_end.map {
			return [it[0],it[1][0]];
			}
	emit:
		versions = version_ch
		paired_end
		single_end
	}


