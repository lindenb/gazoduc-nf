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
include { REGENIE_MAKE_ANNOT  as MAKE_ANNOT        } from '../../../modules/jvarkit/regeniemakeannot/main.nf'
include { flatMapByIndex                           } from '../../../modules/utils/functions.nf'


workflow REGENIE_MAKE_ANNOT {
    take:
        metadata
        scores  // [meta,scores]
        tsv // [meta,tsv] output of regeniefunctional annot  for example
    main:
        versions = Channel.empty()
        multiqc = Channel.empty()

        MAKE_ANNOT(
			scores,
			tsv
			)
		versions = versions.mix(MAKE_ANNOT.out.versions)

		aaf = MAKE_ANNOT.out.aaf
			.flatMap{row->flatMapByIndex(row,1)}
			.map{meta,f->[meta.plus(id: "annot."+ f.name.replaceAll("\\.aaf\\.txt\\.gz\$","")),f]}
			
		
		setfile = MAKE_ANNOT.out.setfile
			.flatMap{row->flatMapByIndex(row,1)}
			.map{meta,f->[meta.plus(id: "annot."+ f.name.replaceAll("\\.setfile\\.txt\\.gz\$","")),f]}
		
		annot = MAKE_ANNOT.out.annot
			.flatMap{row->flatMapByIndex(row,1)}
			.map{meta,f->[meta.plus(id: "annot."+ f.name.replaceAll("\\.annot\\.txt\\.gz\$","")),f]}

		mask = MAKE_ANNOT.out.mask
			.flatMap{row->flatMapByIndex(row,1)}
			.map{meta,f->[meta.plus(id: "annot."+ f.name.replaceAll("\\.mask\\.txt\\.gz\$","")),f]}
		
		annotations = annot.join(setfile).join(aaf).join(mask)



    emit:
        versions
        multiqc
        annotations // [meta,annot,sefile,aaf,mask]
}   
