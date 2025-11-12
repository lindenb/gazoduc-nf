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
include { DOWNLOAD_GTF_OR_GFF3  as DOWNLOAD_GFF3   }  from '../../../modules/gtf/download'
workflow GFF3_INPUT {
take:
    metadata
    fasta
    fai
    dict
main:
    versions = Channel.empty()
    otherwise_ch = Channel.empty()
    gff3_out = Channel.empty()
   
    
    if(metadata.arg_name==null) {
        log.warn("undefined GFF3_INPUT.medata.arg_name")
        }

    if(metadata.require_index==null) {
        log.warn("undefined GFF3_INPUT.medata.require_index")
        }
    if(metadata.download==null) {
        log.warn("undefined GFF3_INPUT.medata.download")
        }
 

    /** DOWNLOAD IF  path is null */
    if(isBlank(metadata.path) && parseBoolean(metadata.download) ) {
        log.warn("--${metadata.arg_name} undefined, trying to download one from fasta...")
        DOWNLOAD_GFF3(
            fasta,
            fai,
            dict
            )
        versions = versions.mix(DOWNLOAD_GFF3.out.versions)
        gff3_out = DOWNLOAD_GFF3.out.gtf//tbi added later
            

        if(!parseBoolean(metadata.require_index))
            {
            gff3_out = gff3_out.map{meta,gff3,tbi->[meta,gff3]}
            }
        }
    else
        {
        if(metadata.path==null ) {
            throw new IllegalArgumentException("GFF3_INPUT  is missing for --${metadata.arg_name}");
            }
       
        else if(metadata.path.endsWith(".gff") || 
                metadata.path.endsWith(".gff3")  ||
                metadata.path.endsWith(".gff.gz")  || 
                metadata.path.endsWith(".gff3.gz") ||
                ( metadata.path.endsWith(".tbi") &&  metadata.path.contains("##idx##")) 
                ) {
                ch1 = Channel.of(metadata.path).map{[gff3:it]}
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
                throw new IllegalArgumentException("Undefined suffix for GFF3: ${metadata.path}");
                }



        // use samtools ##idx## system to split gff3 and index 
        ch1 = ch1.map{
            if(!isBlank(it.tbi)) return it;
            def delim= "##idx##";
            def idx = it.gff3.indexOf(delim);
            if(idx>0) {
                def fa = it.gff3.substring(0,idx);
                def fb = it.gff3.substring(idx+delim.length());
                return it.plus([gff3:fa,tbi:fb]);
                }
            return it;
            }


        ch1 = ch1.map{assertKeyExistsAndNotEmpty(it,"gff3")}
                .map{assertKeyMatchRegex(it,"gff3",".*\\.(gff|gff3|gff\\.gz|gff3\\.gz)")}

       
  

        if(parseBoolean(metadata.require_index?:true)) {
            ch1 = ch1
                .map{
                    if(!it.gff3.endsWith(".gz")) throw new IllegalArgumentException("gff3 fmust have suffix .gz: ${it}");
                    return it;
                    }
                .map{
                    if(isBlank(it.tbi)) {
                        return it.plus(tbi: it.gff3+".tbi");
                        }
                    return it;
                    }
            }    
        
        ch1.count()
            .filter{it!=1}
            .map{
                def msg = "Expected only one GFF for --${metadata.arg_name?:"gff3"} but got N=${it}.";
                log.error(msg);
                throw new IllegalArgumentException(msg);
                }
      

        if(parseBoolean(metadata.require_index?:true)) {
            gff3_out = ch1.map{[
                it.plus([id: removeCommonSuffixes(file(it.gff3).name)]).findAll{k,v->!k.matches("(gff3|tbi)")},
                file(it.gff3),
                file(it.tbi)
                ]}
            otherwise_ch = Channel.of([[id:"no_gff3"],[],[]])
            }
        else
            {
            gff3_out = ch1.map{[
                it.plus([id: removeCommonSuffixes(file(it.gff3).name)]).findAll{k,v->!k.matches("(gff3|tbi)")},
                file(it.gff3)
                ]}
            otherwise_ch = Channel.of([[id:"no_gff3"],[]])
            }

         }

    gff3_out  = gff3_out.first()
emit:
    versions
    gff3 = gff3_out
}