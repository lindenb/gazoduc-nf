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
include { SNPEFF_DOWNLOAD                          } from '../../modules/snpeff/download'
include { SNPEFF_APPLY                             } from '../../modules/snpeff/apply'


workflow SNPEFF {
    take:
        metadata
        snpeff_dictory // meta,snpeff_dictory /* null but not empty */
        vcfs
    main:
        versions = Channel.empty()
        multiqc = Channel.empty()
        SNPEFF_DOWNLOAD(snpeff_dictory)
        versions = versions.mix(SNPEFF_DOWNLOAD.out.versions)

        SNPEFF_APPLY(
            SNPEFF_DOWNLOAD.out.database,
            vcfs
            )
        versions = versions.mix(SNPEFF_APPLY.out.versions)
	    vcfs = 	SNPEFF_APPLY.out.vcf

    emit:
        vcf = SNPEFF_APPLY.out.vcf
        multiqc
        versions
}