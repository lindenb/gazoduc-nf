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
include { isBlank                    }  from '../../../modules/utils/functions.nf'
include { parseBoolean               }  from '../../../modules/utils/functions.nf'
include { assertKeyExistsAndNotEmpty }  from '../../../modules/utils/functions.nf'
include { assertKeyMatchRegex        }  from '../../../modules/utils/functions.nf'
include { removeCommonSuffixes       }  from '../../../modules/utils/functions.nf'
include { htslibSplitIndex           }  from '../../../modules/utils/functions.nf'



workflow VCF_INPUT {
take:
    metadata
main:
    versions = Channel.empty()

    if(metadata.path==null || isBlank(metadata.path)) {
        log.warn("VCF_INPUT : undefined path")
        }

    if(metadata.require_index==null) {
        log.warn("VCF_INPUT : undefined require_index")
        }

    if(isBlank(metadata.arg_name) || isBlank(metadata.arg_name) ) {
        log.warn("VCF_INPUT : undefined arg_name")
        }
    if(metadata.required==null) {
        log.warn("VCF_INPUT : undefined required")
        }
    if(metadata.unique==null) {
        log.warn("VCF_INPUT : undefined unique")
        }

    ch1 = Channel.empty()

    if(metadata.path==null) {
            throw new IllegalArgumentException("VCF_INPUT.meta.path is missing for ${metadata.arg_name?:"option"}");
            }
    else if((metadata.path.endsWith(".tbi")  ||  metadata.path.endsWith(".csi")) && htslibSplitIndex(metadata.path).size()==2) {
            ch1 =  Channel.of(metadata.path).
                    map{htslibSplitIndex(it)}.
                    map{[vcf:it[0],index:it[1]]}
            }
    else if(metadata.path.endsWith(".vcf") || 
            metadata.path.endsWith(".vcf.gz")  ||
            metadata.path.endsWith(".vcf.bgz")  || 
            metadata.path.endsWith(".bcf") ) {
            ch1 = Channel.of(metadata.path).map{[vcf:it]}
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
            throw new IllegalArgumentException("Undefined suffix for VCF: ${metadata.path}");
            }

    ch1 = ch1
            .map{assertKeyExistsAndNotEmpty(it,"vcf")}
            .map{assertKeyMatchRegex(it,"vcf",".*\\.(vcf|vcf\\.gz|bcf)")}



    if(metadata.require_index==null || parseBoolean(metadata.require_index?:true)) {
        ch1 = ch1
            .map{
                if(it.vcf.endsWith(".vcf")) throw new IllegalArgumentException("vcf fmust have suffix .vcf.gz or .bcf: ${it}");
                return it;
                }
            .map{
                if(isBlank(it.index)) {
                    if(it.vcf.endsWith(".vcf.gz")) return it.plus(index: it.vcf+".tbi");
                    return it.plus(index: it.vcf+".csi");
                    }
                return it;
                }
        }



    if(metadata.required==null || parseBoolean(metadata.required?:true)) {
        ch1.count()
            .filter{it==0}
            .map{
                def msg = "Expected at least one VCF for ${metadata.arg_name} but got N=${it}.";
                log.error(msg);
                throw new IllegalArgumentException(msg);
                }
        }

    if(metadata.unique==null || parseBoolean(metadata.unique?:false)) {
        ch1.count()
            .filter{it>1}
            .map{
                def msg = "Expected only one VCF for ${metadata.arg_name} but got N=${it}.";
                log.error(msg);
                throw new IllegalArgumentException(msg);
                }
        }
    
    if(metadata.require_index==null || parseBoolean(metadata.require_index?:true)) {
        ch1 = ch1.map{[
            it.plus([id: removeCommonSuffixes(file(it.vcf).name)]).findAll{k,v->!k.matches("(vcf|index)")},
            file(it.vcf),
            file(it.index)
            ]}
        otherwise_ch = Channel.of([[id:"no_vcf"],[],[]])
        }
    else
        {
        ch1 = ch1.map{[
            it.plus([id: removeCommonSuffixes(file(it.vcf).name)]).findAll{k,v->!k.matches("(vcf|index)")},
            file(it.vcf)
            ]}
        otherwise_ch = Channel.of([[id:"no_vcf"],[]])
        }

    if(metadata.unique!=null || arseBoolean(metadata.unique?:false)) {
        ch1  = ch1.first()
        }

emit:
    versions
    vcf = ch1
}
