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
include {BEDTOOLS_INTERSECT                  } from '../../../modules/bedtools/intersect'
include {BEDTOOLS_SUMMARY                    } from '../../../modules/bedtools/summary'
include {IF_EMPTY                            } from '../../../subworkflows/nf/if_empty'
/**
 * Check user file, if defined, is an ASCII bed, have the correct contigs
 */
workflow PREPARE_USER_BED {
take:
    meta
    fasta
    fai
    dict
    scatter_bed // [meta,bed] genome processed with gatk Scatter
    bed //[meta,userbed] or null or [meta,[]]
main:
    versions = Channel.empty()
    multiqc = Channel.empty()

        CHECK_BED(
            fasta,
            fai,
            dict,
            bed
            )
        versions = versions.mix(CHECK_BED.out.versions)
        bed = CHECK_BED.out.bed

        if(meta.with_intersect!=null && meta.with_intersect!=false) {
            BEDTOOLS_INTERSECT(
                    PREPARE_REFERENCE.out.fai,
                    PREPARE_REFERENCE.out.scatter_bed.combine(bed)
                    )
            versions = versions.mix(BEDTOOLS_INTERSECT.out.versions)
            bed  = BEDTOOLS_INTERSECT.out.bed.first()
            }



    	bed = IF_EMPTY(bed, scatter_bed)

    BEDTOOLS_SUMMARY(fai,bed)
    versions = versions.mix(BEDTOOLS_SUMMARY.out.versions)
    multiqc = multiqc.mix( BEDTOOLS_SUMMARY.out.summary )
emit:
    versions
    bed
    multiqc
}

process CHECK_BED {
label "process_short"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
tag "${meta.id?:bam.name}"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
	tuple val(meta ),path(bed)
output:
	tuple val(meta ),path("*.bed"),emit:bed
	path("versions.yml"),emit:versions
script:
	def convert = task.ext.convert?:""
    def prefix = task.ext.prefix?:"${meta.id}.checked"
    def with_merge = (task.ext.with_merge?:false) as boolean
"""
	hostname 1>&2
    set -x
	mkdir -p TMP
    
    ${bed.name.endsWith(".gz") || bed.name.endsWith(".bgz")?"gunzip -c":"cat"} "${bed}" > TMP/jeter.bed
	
    
    file  TMP/jeter.bed | cut -d ':' -f2- | tee TMP/type.txt 1>&2
    # the following test wills FAIL if type doesn't contains "ASCII text"
    grep -F 'ASCII text' TMP/type.txt 1>&2

    # the following test wills FAIL if type contains CRLF, LF line terminators
    # see https://en.wikipedia.org/wiki/Newline#Issues_with_different_newline_formats
    grep -L -F 'CRLF, LF' TMP/type.txt > TMP/jeter.err || true
    test -s TMP/jeter.err
    
    # cleanup headers
    grep -v '^(#|browser|track)' TMP/jeter.bed > TMP/jeter2.bed
    mv TMP/jeter2.bed TMP/jeter.bed

    # the following test wills FAIL there is a line without 3 columns
    awk -F '\t' '(\$NF <3)' TMP/jeter.bed > TMP/jeter.err
    head TMP/jeter.err 1>&2
    test ! -s TMP/jeter.err

    # the following test wills FAIL if start or end != integer
    cut -f2,3 TMP/jeter.bed | tr "\t" "\\n" |\\
        grep -v  -x -E  '[0-9]+' > TMP/jeter.err || true
    head TMP/jeter.err | awk '{printf("NOT AN INTEGER: (%s)\\n",\$1);}' 1>&2
    test ! -s TMP/jeter.err

    # the following test wills FAIL there is a line with start>=end
    awk -F '\t' '(\$2 >= \$3)' TMP/jeter.bed > TMP/jeter.err
    head TMP/jeter.err 1>&2
    test ! -s TMP/jeter.err

    # rename chromosomes
	jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP bedrenamechr \\
            -R ${fasta} \\
            -c 1  \\
            ${convert.isEmpty()?"":"--convert ${convert}"} \\
            TMP/jeter.bed > TMP/jeter2.bed
    mv TMP/jeter2.bed TMP/jeter.bed
	
    # sort bed
    sort -T TMP -t '\t' -k1,1 -k2,2n TMP/jeter.bed  > TMP/jeter2.bed
    mv TMP/jeter2.bed TMP/jeter.bed

    # merge if enabled
    if ${with_merge}
    then
        bedtools merge -i TMP/jeter.bed |\\
            sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter2.bed
        mv TMP/jeter2.bed TMP/jeter.bed
    fi

    # the following test wills FAIL if any END in bed is longer than chromosome
    join -t '\t' -1 1 -2 1 -o '1.1,1.2,2.2' \\
        <(cut -f1,3 TMP/jeter.bed | sort -T TMP -t '\t' -k1,1) \\
        <(cut -f1,2 "${fai}" | sort -T TMP -t '\t' -k1,1) |\\
        awk -F '\t' '(int(\$2) > int(\$3))' > TMP/jeter.err
    head TMP/jeter.err 1>&2
    test ! -s TMP/jeter.err


    mv TMP/jeter.bed ${prefix}.bed

cat << END_VERSIONS > versions.yml
"${task.process}":
	jvarkit: todo
END_VERSIONS
"""

stub:
def prefix = task.ext.prefix?:"${meta.id}.checked"
"""
touch versions.yml cp "${bed}" "${prefix}.bed"
"""
}
