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
include { isBlank                                  }  from '../../../modules/utils/functions.nf'
include { parseBoolean                             }  from '../../../modules/utils/functions.nf'
include { assertKeyExistsAndNotEmpty               }  from '../../../modules/utils/functions.nf'
include { assertKeyMatchRegex                      }  from '../../../modules/utils/functions.nf'
include { removeCommonSuffixes                     }  from '../../../modules/utils/functions.nf'
include { DOWNLOAD_GTF_OR_GFF3   as DOWNLOAD_GTF   }  from '../../../modules/gtf/download'
include { htslibSplitIndex                         }  from '../../../modules/utils/functions.nf'

workflow GTF_INPUT {
take:
    metadata
    dict
main:
    versions = Channel.empty()
    otherwise_ch = Channel.empty()
    gtf_out = Channel.empty()
   
    
    if(metadata.arg_name==null) {
        log.warn("undefined GTF_INPUT.medata.arg_name")
        }

    if(metadata.require_index==null) {
        log.warn("undefined GTF_INPUT.medata.require_index")
        }
    if(metadata.download==null) {
        log.warn("undefined GTF_INPUT.medata.download")
        }
 

    /** DOWNLOAD IF  path is null */
    if(isBlank(metadata.path) && parseBoolean(metadata.download) ) {
        log.warn("--${metadata.arg_name} undefined, trying to download one from dict...")
        DOWNLOAD_GTF(dict)
        versions = versions.mix(DOWNLOAD_GTF.out.versions)
        gtf_out = DOWNLOAD_GTF.out.gtf//tbi added later
            

        if(!parseBoolean(metadata.require_index))
            {
            gtf_out = gtf_out.map{meta,gtf,tbi->[meta,gtf]}
            }

        }
    else
        {
        if(metadata.path==null ) {
            throw new IllegalArgumentException("GTF_INPUT  is missing for --${metadata.arg_name}");
            }
       
        else if(metadata.path.endsWith(".gtf") || 
                metadata.path.endsWith(".gtf.gz") ||
                ( metadata.path.endsWith(".tbi") &&  htslibSplitIndex(metadata.path).size()==2) 
                ) {
                ch1 = Channel.of(metadata.path).map{[gtf:it]}
                }
        else if(metadata.path.endsWith(".list") || metadata.path.endsWith(".txt")) {
                ch1 = Channel.fromPath(metadata.path)
                    .splitText()
                    .map{it.trim()}
                    .filter{!isBlank(it)}
                }
        else if(metadata.path.endsWith(".json")) {
                ch1 = Channel.fromPath(metadata.path)
                    .splitJson()
                }
        else if(metadata.path.endsWith(".csv")) {
                ch1 = Channel.fromPath(metadata.path)
                    .splitCsv(header:true,sep:',')
                }
        else if(metadata.path.endsWith(".tsv")) {
                ch1 = Channel.fromPath(metadata.path)
                    .splitCsv(header:true,sep:'\t')
                }
        else
                {
                throw new IllegalArgumentException("Undefined suffix for GTF: ${metadata.path}");
                }



        // use samtools ##idx## system to split gtf and index 
        ch1 = ch1.map{
             if(!isBlank(it.tbi)) return it;
            def a =  htslibSplitIndex(it.gtf);
            if(a.size()==2) {
                 return it.plus([gtf:a[0],tbi:a[1]]);
                }
            return it;
            }


        ch1 = ch1.map{assertKeyExistsAndNotEmpty(it,"gtf")}
                .map{assertKeyMatchRegex(it,"gtf",".*\\.(gtf|gtf\\.gz)")}

       
  

        if(parseBoolean(metadata.require_index?:true)) {
            ch1 = ch1
                .map{
                    if(!it.gtf.endsWith(".gz")) throw new IllegalArgumentException("gtf fmust have suffix .gz: ${it}");
                    return it;
                    }
                .map{
                    if(isBlank(it.tbi)) {
                        return it.plus(tbi: it.gtf+".tbi");
                        }
                    return it;
                    }
            }    
        
        ch1.count()
            .filter{it!=1}
            .map{
                def msg = "Expected only one GTF for --${metadata.arg_name?:"gtf"} but got N=${it}.";
                log.error(msg);
                throw new IllegalArgumentException(msg);
                }
      

        if(parseBoolean(metadata.require_index?:true)) {
            gtf_out = ch1.map{[
                it.plus([id: removeCommonSuffixes(file(it.gtf).name)]).findAll{k,v->!k.matches("(gtf|tbi)")},
                file(it.gtf),
                file(it.tbi)
                ]}
            otherwise_ch = Channel.of([[id:"no_gtf"],[],[]])
            }
        else
            {
            gtf_out = ch1.map{[
                it.plus([id: removeCommonSuffixes(file(it.gtf).name)]).findAll{k,v->!k.matches("(gtf|tbi)")},
                file(it.gtf)
                ]}
            otherwise_ch = Channel.of([[id:"no_gtf"],[]])
            }

         }

    gtf_out  = gtf_out.first()
emit:
    versions
    gtf = gtf_out
}
