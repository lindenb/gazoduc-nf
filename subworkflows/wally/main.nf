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
include { DGV_DOWNLOAD                     } from '../../modules/dgv/download'
include { verify                           } from '../../modules/utils/functions.nf'
include { isBlank                          } from '../../modules/utils/functions.nf'
include { flatMapByIndex                   } from '../../modules/utils/functions.nf'
include { DOWNLOAD_GNOMAD_SV               } from '../../modules/gnomad_sv/download.vcf'
include { BCFTOOLS_QUERY                   } from '../../modules/bcftools/query'
include { WALLY_REGION                     } from '../../modules/wally/region'
include { ZIP as ZIP_BY_INTERVAL           } from '../../modules/utils/zip'
include { ZIP as ZIP_ALL                   } from '../../modules/utils/zip'

workflow WALLY {
    take:
        metadata
        fasta
        fai
        dict
        bams // Channel[meta,bam,bai] //meta contains contig/start/end
    main:
        versions = Channel.empty()
        multiqc = Channel.empty()

        DGV_DOWNLOAD(dict)
        versions = versions.mix(DGV_DOWNLOAD.out.versions)

        DOWNLOAD_GNOMAD_SV(dict)
        versions = versions.mix(DOWNLOAD_GNOMAD_SV.out.versions)

        BCFTOOLS_QUERY(DOWNLOAD_GNOMAD_SV.out.vcf.map{meta,vcf,_tbi->[meta,vcf]})
        versions = versions.mix(BCFTOOLS_QUERY.out.versions)
        
        MERGE_KNOWN(
            DGV_DOWNLOAD.out.bed
                .map{meta,bed,tbi->bed}
                .mix(BCFTOOLS_QUERY.out.output.map{meta,bed->bed})
                .collect()
                .map{files->[[id:"known"],files.sort()]}
            )
        versions = versions.mix(MERGE_KNOWN.out.versions)

        WALLY_REGION(
            fasta,
            fai,
            MERGE_KNOWN.out.bed,
            bams.map{meta,bam,bai->[meta,bam,bai,[]]}
            )
        versions = versions.mix(WALLY_REGION.out.versions)

        png_ch =WALLY_REGION.out.png
                .flatMap{row->flatMapByIndex(row,1)}
                .map{meta,png->[[id:"${meta.contig}_${meta.start}_${meta.end}",contig:meta.contig,start:meta.start,end:meta.end],png]}
        html_ch =WALLY_REGION.out.html
                .flatMap{row->flatMapByIndex(row,1)}
                .map{meta,html->[[id:"${meta.contig}_${meta.start}_${meta.end}",contig:meta.contig,start:meta.start,end:meta.end],html]}
                


        HTML_PAGE(html_ch.groupTuple())
        versions = versions.mix(HTML_PAGE.out.versions)

        ZIP_BY_INTERVAL(
            png_ch.mix(HTML_PAGE.out.html).groupTuple()
            )
        versions = versions.mix(ZIP_BY_INTERVAL.out.versions)

        ZIP_ALL(
            ZIP_BY_INTERVAL.out.zip
                .map{meta,z->z}
                .collect()
                .map{files->[metadata,files.sort()]}
            )
        versions = versions.mix(ZIP_BY_INTERVAL.out.versions)

    emit:
        versions
        multiqc 
}


process MERGE_KNOWN {
label "process_single"
tag "${meta.id?:""}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
        tuple val(meta),path("INPUT/*")
output:
        tuple val(meta),path("*.bed.gz"),path("*.bed.gz.tbi"),emit:bed
        path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"${meta.id}"
"""
hostname 1>&2
set -o pipefail

gunzip -c INPUT/*.gz |\
        grep -v "^#" |\
        cut -f 1-4 |\
        awk -F '\t' '(\$2 != \$3)' |\
        sort -t '\t' -T . -k1,1 -k2,2n -k3,3n --unique |\
        uniq > ${prefix}.bed

bgzip ${prefix}.bed
tabix -p bed ${prefix}.bed.gz
touch versions.yml
"""

stub:
    def prefix = task.ext.prefix?:"${meta.id}"
"""
touch versions.yml  ${prefix}.bed  ${prefix}.bed.gz
"""
}

process HTML_PAGE {
label "process_single"
tag "${meta.id}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
        tuple val(meta),path(htmls)
output:
        tuple val(meta),path("*.html"),emit:html
        path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"${meta.id}"
"""
mkdir -p TMP

cat << 'EOF' > TMP/jeter.html
<html>
  <head>
     <title>${meta.contig}:${meta.start}-${meta.end}</title>
   </head>
<body>
<div>
<h1>${meta.contig}:${meta.start}-${meta.end}</h1>
<div>
EOF

cat ${htmls.join(" ")} >> TMP/jeter.html

cat << EOF >> TMP/jeter.html
</div></div></body></html>
EOF

mv TMP/jeter.html  ${prefix}.html

touch versions.yml
"""
stub:
"""
touch versions.yml ${meta.id}.html
"""
}
