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
include { ANNOTSV_INSTALL           } from '../../modules/annotsv/install'
include { ANNOTSV_APPLY             } from '../../modules/annotsv/apply'

workflow ANNOTSV {
    take:
        metadata
        fasta
        fai
        dict
        vcfs
    main:
        versions = Channel.empty()
        multiqc = Channel.empty()
        
        
        ANNOTSV_INSTALL(
            dict
            .filter{meta,_dict->meta.ucsc_name!=null && meta.ucsc_name.matches("hg(19|38)")}
            .map{metadata}
            )
        versions = versions.mix(ANNOTSV_INSTALL.out.versions)
        
        ANNOTSV_APPLY(
            ANNOTSV_INSTALL.out.annotsv_db,
            dict,
            vcfs
            )
        versions = versions.mix(ANNOTSV_APPLY.out.versions)

    emit:
        versions
        multiqc
}