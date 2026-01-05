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

include { GRAPHTYPER as GTYPER                 }  from '../../../modules/graphtyper/genotype'
include { BCFTOOLS_CONCAT                      }  from '../../../modules/bcftools/concat3'
include { verify                               }  from '../../../modules/utils/functions.nf'




workflow GRAPHTYPER {
take:
    metadata
    fasta
    fai
    dict
    mosdepth_summary //[meta, []]
    beds // meta,bed
    bams // meta,bam,bai
main:
    versions = Channel.empty()
    multiqc = Channel.empty()

    // check uniq bed
    beds.map{meta,bed->[meta.id,bed]}
        .groupTuple()
        .map{
            verify(it[1].size()==1,"GRAPHTYPER : bed id should be unique ${it}.");
            return -1;
            }
   
    // create a key as pivot to join mosdepth and bams
    ch1 = mosdepth_summary.map{meta,f->[meta.id,meta,f]}
    ch2 = bams.map{meta,bam,bai->[meta.id, meta,bam,bai]}
    
    ch3 = ch1.join(ch2).//sample, meta1,summary,meta2,bam,bai
        map{_meta_id,meta1,summary,meta2,bam,bai->[meta2,bam,bai,summary]} // meta, bam,bai, summary

    COVERAGE_DIVIDE_READLENGTH(fasta,fai,ch3)
    versions = versions.mix(COVERAGE_DIVIDE_READLENGTH.out.versions)

    merged = COVERAGE_DIVIDE_READLENGTH.out.output
        .map{meta,f->f}
        .collectFile(name: 'lengths.tsv')
        .map{[[id:"gtyper"],it]}
        .first()

    
    bam_bed_ch = beds
        .combine(bams.map{_meta,bam,bai->[bam,bai]}.flatMap().collect().map{it.sort()}.toSortedList())
        .map{meta,bed,bams->[meta,bams,bed]}//meta, BAMS, bed
        

    GTYPER(
        fasta,
        fai,
        merged,
        bam_bed_ch
        )
    versions = versions.mix(GTYPER.out.versions)

    BCFTOOLS_CONCAT(
        GTYPER.out.vcf
            .map{_meta,vcf,tbi->[vcf,tbi]}//gvcf,tbi
             .flatMap()
             .collect()
             .map{[[id:"graphtyper"],it.sort()]},
        )
    versions = versions.mix(BCFTOOLS_CONCAT.out.versions)
    

emit:
    versions
    vcf = BCFTOOLS_CONCAT.out.vcf
    multiqc
}




process COVERAGE_DIVIDE_READLENGTH {
tag "${meta.id?:""}"
label "process_single"
label "array100"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
	tuple val(meta),path(bam),path(bai),path(mosdepth_summary)
output:
	tuple val(meta),path("*.tsv"),emit:output
	path("versions.yml"),emit:versions
script:
	def n_reads = task.ext.n_reads?:1000
    def prefix = task.ext.prefix?:"${meta.id}"
"""
hostname 1>&2
mkdir -p TMP

COV=\$(awk -F '\t' '(\$1=="total_region") {print \$4}' "${mosdepth_summary}") 
test ! -z "\${COV}"

set +o pipefail

samtools view -F 3844 -T "${fasta}" "${bam}" | \\
        head -n "${n_reads}" |\\
		awk -F '\t' -vCOV=\${COV} 'BEGIN{T=0.0;N=0;} {N++;T+=length(\$10)} END{printf("${meta.id}\t%f\\n",COV/(N==0?100:(T/N)));}' > TMP/jeter.tsv

test -s TMP/jeter.tsv

mv TMP/jeter.tsv ${prefix}.tsv

cat << EOF > versions.yml
"${task.process}"
    "samtools": "todo"
EOF
"""
}
