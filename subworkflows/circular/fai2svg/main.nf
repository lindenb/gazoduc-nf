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
include { FAI2BED              } from '../../../modules/samtools/fai2bed'
include { BED_TO_XML           } from '../../../modules/jvarkit/bed2xml'
include { XSLTPROC             } from '../../../modules/xsltproc'

workflow FAI_TO_SVG {
take:
    metadata
    fai
    dict
main:
    versions = Channel.empty()
    multiqc  = Channel.empty()

    FAI2BED(fai)
    versions = versions.mix(FAI2BED.out.versions)

    
    BED_TO_XML(
        dict,
        FAI2BED.out.bed
        )
     versions = versions.mix(BED_TO_XML.out.versions)
   
    XSLTPROC(
        [[id:"no_stylesheet"],"${projectDir}/src/xsl/fai2svg.xsl"],
        BED_TO_XML.out.xml
        )
    versions = versions.mix(XSLTPROC.out.versions)
emit:
    versions
    multiqc
    dict = BED_TO_XML.out.dict
    svg = XSLTPROC.out.xml
}