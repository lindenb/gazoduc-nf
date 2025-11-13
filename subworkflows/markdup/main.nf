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

include { MARK_DUPLICATES as GAKT_MARK_DUPLICATES     } from '../../modules/gatk/markduplicates'
include { SAMTOOLS_MERGE                              } from '../../modules/samtools/merge'
include { SAMTOOLS_MARKDUP                            } from '../../modules/samtools/markdup'
include { SAMBAMBA_MARKDUP                            } from '../../modules/sambamba/markdup'
include { SAMTOOLS_COLLATE                            } from '../../modules/samtools/collate'
include { SAMTOOLS_FIXMATE                            } from '../../modules/samtools/fixmate'
include { SAMTOOLS_SORT                               } from '../../modules/samtools/sort'
include { PB_MARKDUP                                  } from '../../modules/parabricks/markdup'


workflow MARK_DUPLICATES {
	take:
		meta
        fasta
        fai
        dict
		bams // [meta,bam] or [meta,bam,bai]
	main:
		versions = Channel.empty()
		out_bams = Channel.empty()
		multiqc = Channel.empty()

		if(meta.markdup_enabled==null) {
            log.warn("MARKDUP undefined meta.markdup_enabled")
            }
		
	
		if(meta.markdup_enabled==null || meta.markdup_enabled==true) {
            if(meta.markdup_method==null) {
                log.warn("MARKDUP undefined meta.markdup_method")
                }
            

        if(meta.markdup_method==null || meta.markdup_method.matches("(markduplicates|gatk|picard|pb|parabricks)")) {
            GAKT_MARK_DUPLICATES(
                fasta,
                fai,
                dict,
                bams
                    .map{it.size()==3?it: [it[0],it[1],null]} // bai provided ?
                    .map{meta,bam,bai->[meta.id,meta,bam]}
                    .groupTuple()
                    .map{id,metas,bam_files->[metas[0],bam_files.sort()]}
                )
            versions = versions.mix(GAKT_MARK_DUPLICATES.out.versions)
            out_bams = GAKT_MARK_DUPLICATES.out.bam
            multiqc = multiqc.mix(GAKT_MARK_DUPLICATES.out.metrics)
            }
        else if(meta.markdup_method.equals("sambamba") || meta.markdup_method.equals("samtools")) {
            
            /* branch so bams that have been splitted  */
            ch1 =   bams
                .map{it.size()==3?it: [it[0],it[1],null]} // bai provided
                .map{meta,bam,bai->[meta.id,meta,bam]}
                .groupTuple()
                .map{id,metas,bam_files->[metas[0],bam_files.sort()]}
                .branch{
                    needmerge: it[1].size() > 1
                    one: true
                    }
            
            SAMTOOLS_MERGE(
                fasta,
                fai,
                ch1.needmerge.map{meta,bam_files->[meta,bam_files.sort()]}
                )
            versions = versions.mix(SAMTOOLS_MERGE.out.versions)
            

            merged_ch =  SAMTOOLS_MERGE.out.bam
                    .map{meta,bam,bai->[meta,bam]}
                    .mix(
                        ch1.one.map{meta,bam_files->[meta,bam_files[0]]}
                        )

            if(meta.markdup_method.equals("sambamba")) {
                SAMBAMBA_MARKDUP(merged_ch)
                versions = versions.mix(SAMBAMBA_MARKDUP.out.versions)
                out_bams = SAMBAMBA_MARKDUP.out.bam
                }
            else
                {
                SAMTOOLS_COLLATE(fasta,fai,merged_ch)
                versions = versions.mix(SAMTOOLS_COLLATE.out.versions)

                SAMTOOLS_FIXMATE(fasta,fai,SAMTOOLS_COLLATE.out.bam)
                versions = versions.mix(SAMTOOLS_FIXMATE.out.versions)

                SAMTOOLS_SORT(fasta,fai,SAMTOOLS_FIXMATE.out.bam)
                versions = versions.mix(SAMTOOLS_SORT.out.versions)


                SAMTOOLS_MARKDUP(fasta,fai,SAMTOOLS_SORT.out.bam)
                versions = versions.mix(SAMTOOLS_MARKDUP.out.versions)
                out_bams = SAMTOOLS_MARKDUP.out.bam
                multiqc = multiqc.mix(SAMTOOLS_MARKDUP.out.json)
                }
            }
        else if(meta.markdup_method.matches("(pb|parabricks)"))  {
            PB_MARKDUP(
                fasta,
                fai,
                merged_ch
                )
            versions = versions.mix(PB_MARKDUP.out.versions)
            multiqc = multiqc.mix(PB_MARKDUP.out.metrics)
            }
        else
            {
            throw new IllegalArgumentException("Boum MAP_BWA undefined meta.markdup_method: ${meta.markdup_method}");
            }
        }
    else
        {
        out_bams = bams
        }

emit:
    versions
    bams = out_bams
    multiqc
}
