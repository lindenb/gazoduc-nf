/*

Copyright (c) 2026 Pierre Lindenbaum

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
include {assertKeyExistsAndNotEmpty          } from '../../../modules/utils/functions.nf'
include {assertKeyMatchRegex                 } from '../../../modules/utils/functions.nf'
include {isBlank                             } from '../../../modules/utils/functions.nf'
include {SAMTOOLS_SAMPLES                    } from '../../../modules/samtools/samples'
include {GROUP_BY_TRIOS                      } from '../../../subworkflows/pedigree/groupbytrios'




/** decode meta to bams, allow multiple fasta references */
workflow META_TO_BAMS {
take:
	metadata
    fasta // [meta,fasta] MAIN REFERENCE FASTA
    fai   // [meta,fai  ] MAIN REFERENCE FAI
    dict  // [meta,dict  ] MAIN REFERENCE DICT
	metas
main:
	versions = Channel.empty()

    //first create a [META,BAM,BAI]
    bams = metas
        .map{assertKeyMatchRegex(it,"bam","^\\S+\\.(bam|cram)\$")}
        .map{
            if(!isBlank(it.bai)) return it;
            if(it.bam.endsWith(".cram")) return it.plus(bai : it.bam+".crai");
            return it.plus(bai:it.bam+".bai");
        	}
		.map{assertKeyMatchRegex(it,"bai","^\\S+\\.(bai|crai)\$")}
    	.map{[
            it.findAll{k,v->!k.matches("(bam|bai)")},
            file(it.bam),
            file(it.bai)
            ]}

    // find those having a defined FASTA
    ch1 = bams.branch{v->
        without_ref: isBlank(v[0].fasta)
        with_ref   : true
        }

    // main ref linearized
    main_ref_ch = fasta
        .combine(fai)
        .combine(dict)
        .map{meta1,fasta,_meta2,fai,_meta3,dict->[meta1,fasta,fai,dict]}
        .first()

    // add default fasta to bams missing a ref
    without_ref = ch1.without_ref
        .combine(main_ref_ch)
        .map{meta1,bam,bai,_meta2,fasta,fai,dict->[meta1,bam,bai,fasta,fai,dict]}
    
    // add  meta.fasta to bams having a ref
    with_ref = ch1.with_ref
        .map{meta,bam,bai->[meta,bam,bai,file(meta.fasta)]}
        .map{meta,bam,bai,fasta->
            def parent = fasta.getParent();
            if(parent==null) throw new IllegalArgumentException("cannot get parent of fasta ${fasta}");
            def fai = fasta.getParent().resolve(fasta.name +".fai");
            return [meta,bam,bai,fasta,fai];
            }
        .map{meta,bam,bai,fasta,fai->
            def parent = fasta.getParent();
            if(parent==null) throw new IllegalArgumentException("cannot get parent of fasta ${fasta}");
            def dict= fasta.getParent().resolve(fasta.getBaseName()+".dict");
            return [
                meta.findAll{k,v->!k.matches("(fasta|fai|dict)")},
                bam,
                bai,
                fasta,
                fai,
                dict
                ];
            }
    
    bams = without_ref.mix(with_ref);


    /* find samples without sample */
    ch2 = bams.branch{v->
        no_sample: isBlank(v[0].sample) 
        with_sample: true
        }

    
    SAMTOOLS_SAMPLES(
        main_ref_ch
            .flatMap{meta,fasta,fai,dict->[fasta,fai]}
            .collect()
            .map{[[id:"ref"],it.sort()]},
        ch2.no_sample
            .map{meta,bam,bai,fasta,fai,dict->bam}
            .collect()
            .map{[[id:"bams"],it.sort()]}
        )
    
    versions = versions.mix(SAMTOOLS_SAMPLES.out.versions)

        
    ch3 = SAMTOOLS_SAMPLES.out.samplesheet
        .map{it[1]}
        .splitCsv(sep:'\t',header:false)
        .map{if(isBlank(it[0])) throw new IllegalArgumentException("No sample name for ${it}"); return it;}
        .map{row->[file(row[1]).toRealPath().toString(),[sample:row[0]]]}
    

    fixed_sample_ch= ch2.no_sample
        .map{meta,bam,bai,fasta,fai,dict->[bam.toRealPath().toString(),meta,bam,bai,fasta,fai,dict]}
        .join(ch3, failOnDuplicate:true, failOnMismatch:true)
        .map{bam_path,meta1,bam,bai,fasta,fai,dict,meta2->[meta2.plus(meta1),bam,bai,fasta,fai,dict]}
        
        
    // merge fixed sample, append meta.id
    bams_out = ch2.with_sample
        .mix(fixed_sample_ch)
        .map{meta,bam,bai,fasta,fai,dict->[
            (isBlank(meta.id)?meta.plus(id:meta.sample):meta),
            bam,
            bai,
            fasta,
            fai,
            dict
            ]}
    


    if(metadata.with_test_unique_id==null || metadata.with_test_unique_id==true) {
        bams_out.map{it->it[0].id}.unique().count().
            combine(bams_out.map{it->it[0].id}.count())
            .filter{c1,c2->c1!=c2}
            .map{throw new IllegalArgumentException("input contains same duplicate samples");}
        }

    if(metadata.enable_group_by_trios==null || metadata.enable_group_by_trios==true) {
        GROUP_BY_TRIOS(metadata,bams_out)
        bams_out = GROUP_BY_TRIOS.out.output
        versions = versions.mix(GROUP_BY_TRIOS.out.versions)
        }

    all_references =  bams_out
        .flatMap{meta,bam,bai,fasta,fai,dict->[fasta,fai,dict]}
        .mix(fasta.map{it[1]})
        .mix(fai.map{it[1]})
        .mix(dict.map{it[1]})
        .filter{fn->fn.exists()} // when running in stub mode...
        .map{fn->[fn.name,fn.toRealPath()]} // group files by names. prevent file collisiton; FAI might have same name because PREPARE_ONE_REFERENCE.out.fai
	.groupTuple()
	.map{
		def L0 = it[1].collect{f->f.toRealPath().toString()}.unique()
		if(L0.size()!=1) throw new IllegalArgumentException("META2BAMS2: your using several FASTA/FAI/DICT sequences with the same name: ${it[1]}.");
		return it;
		}
        .map{name,fns->fns.sort()[0]}
        .unique()
        .collect()
        .map{[metadata,it]} 

emit:
	versions
	all_references
	bams = bams_out // meta,bam,bai,fasta,fai,dict

}
