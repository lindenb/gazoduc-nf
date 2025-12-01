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
include {isBlank } from '../../modules/utils/functions.nf'
include {verify  } from '../../modules/utils/functions.nf'

workflow CNVNATOR {
take:
    metadata
    fasta
    fai
    dict
    bams // tuple(meta,bam,bai)
main:
    versions  = Channel.empty()
    multiqc = Channel.empty()
    ONE_FILE_PER_CONTIG(fasta,fai)
    versions = versions.mix(ONE_FILE_PER_CONTIG.out.versions)

    each_contig_ch = ONE_FILE_PER_CONTIG.out.chroms_list.
        map{meta,f->f}.
        splitText().
        map{it.trim()}

    
    CALL_CNV(
        fasta,
        fai,
        ONE_FILE_PER_CONTIG.out.chromdir,
        each_contig_ch.combine(bams)
        )
    versions = versions.mix(CALL_CNV.out.versions)

    CNVNATOR_CONCAT_VCFS(
        CALL_CNV.out.vcf
            .groupTuple()
            .map{meta,vcfs,idxs->[meta,vcfs.plus(idxs).flatten().sort()]}
        )
    versions = versions.mix(CNVNATOR_CONCAT_VCFS.out.versions)

    CNVNATOR_CONCAT_BEDS(
        CALL_CNV.out.bed.
            map{it[1]}.
            collect().
            map{[[id:"cnvnator"],it]}
        )
    versions = versions.mix(CNVNATOR_CONCAT_BEDS.out.versions)
emit:
    versions
    multiqc

}



process ONE_FILE_PER_CONTIG {
tag "${fasta.name}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
label "process_single"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
output:
    tuple val(meta1),path("CHROMS"),emit:chromdir
    tuple val(meta1),path("CHROMS/chroms.txt"),emit:chroms_list
    path("versions.yml"),emit:versions
script:
    def regex = task.ext.regex?:"^(chr)?[0-9XY]+\$"
"""
mkdir -p TMP
cut -f 1 "${fai}" | grep -E '${regex}' | while read C
do
        samtools faidx "${fasta}" "\${C}" | tr "\t" " " | cut -f 1 -d ' ' > "TMP/\${C}.fa"
        samtools faidx "TMP/\${C}.fa"
        echo "\${C}" >> TMP/chroms.txt
done

mv TMP CHROMS

cat << END_VERSIONS > versions.yml
${task.process}:
    samtools: \$(samtools version | awk '(NR==1)  {print \$NF}')
END_VERSIONS
"""

stub:
"""
mkdir -p CHROMS
touch versions.yml
cut -f1 '${fai}' > CHROMS/chroms.txt
"""
}




process CALL_CNV {
tag "${meta.id} ${contig} ${bam.name}"
label "process_single"
label "array100"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/cnvnator.yml"

input:
        tuple val(meta1),path(fasta)
        tuple val(meta2),path(fai)
        tuple val(meta4),path(contigdir) 
        tuple val(contig),val(meta),path(bam),path(bai)
output:
        tuple val(meta),path("*.bcf"),path("*.bcf.csi"),emit:vcf
        tuple val(meta),path("*.bed.gz"),emit:bed
        path("versions.yml"),emit:versions
script:
        verify(!isBlank(meta1.ucsc_name),"undefined ucsc_name in ${meta1}?");
        def bin = task.ext.bin?:5000
        def prefix = (task.ext.prefix?:meta.id)+"."+contig+"."+bin
"""
mkdir -p TMP/REF
set -x
set -o pipefail

if ${bam.name.endsWith(".cram")} ; then
        samtools view  --threads ${task.cpus} -O BAM -o TMP/jeter.bam -F 3844 --reference "${fasta}" "${bam}" "${contig}"
        samtools index --threads ${task.cpus} TMP/jeter.bam
fi

        cnvnator \\
                -root "TMP/tmp.root" \\
                -chrom ${contig} \\
                -genome '${meta1.usc_name}' \\
                -tree ${bam.name.endsWith(".cram")?"TMP/jeter.bam":"\"${bam}\""} 1>&2


        # generate histogram
        cnvnator \\
                -root "TMP/tmp.root" \\
                -d "${contigdir}" \\
                -his ${bin} 1>&2
 
        # calculate statistics
        cnvnator \
                -root "TMP/tmp.root" \\
                -d "${contigdir}" \\
                -stat ${bin}  1>&2

        # partition
        cnvnator \
                -root "TMP/tmp.root" \\
                -d "${contigdir}" \\
                -partition ${bin} 1>&2

        # call CNVs
        cnvnator \
                -root "TMP/tmp.root" \\
                -d "${contigdir}" \\
                -call ${bin} > "TMP/output.tsv" 


        cnvnator2VCF.pl \\
                -prefix "${bin}" \
                -reference  '${meta1.ucsc_name}' \\
                "TMP/output.tsv" \
                "${contigdir}" > TMP/jeter.vcf

        # rename header
        bcftools query -l TMP/jeter.vcf |  awk '{printf("%s\\t${meta.id}\\n",\$1);}' > TMP/jeter.tsv
        bcftools reheader --fai "${fai}" --samples TMP/jeter.tsv TMP/jeter.vcf > TMP/jeter2.vcf
        mv TMP/jeter2.vcf TMP/jeter.vcf

        #sort
        bcftools sort  -o "TMP/jeter.bcf" -O b -T TMP TMP/jeter.vcf
        bcftools index TMP/jeter.bcf

        head TMP/jeter.tsv
        awk -F '\t' '{split(\$2,a,/[:\\-]/);printf("%s\t%d\t%s\t%s\t%s\\n",a[1],int(a[2])-1,a[3],"${meta.id}",\$0);}' TMP/output.tsv |\
                grep -v 'Number of free parameters' |\\
                sort -T TMP -t '\t' -k1,1 -k2,2n |\\
                gzip  > TMP/jeter.bed.gz


        mv TMP/jeter.bed.gz ${prefix}.bed.gz
        mv TMP/jeter.bcf ${prefix}.bcf
        mv TMP/jeter.bcf.csi ${prefix}.bcf.csi


cat << END_VERSIONS > versions.yml
${task.process}:
    samtools: \$(samtools version | awk '(NR==1)  {print \$NF}')
    bcftools: \$(bcftools version | awk '(NR==1)  {print \$NF}')
    cnvnator: \$(cnvnator 2>&1 | sed -n '3p' | sed 's/CNVnator v//') (bin:${bin})
END_VERSIONS
"""

stub:
        def bin = task.ext.bin?:5000
        def prefix = (task.ext.prefix?:meta.id)+"."+contig+"."+bin
"""
touch versions.yml ${prefix}.bed.gz ${prefix}.bcf ${prefix}.bcf.csi
"""
}

process CNVNATOR_CONCAT_VCFS {
tag "${meta.id}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
    tuple val(meta),path("VCFS/*")
output:
    tuple val(meta),path("*.bcf"),path("*.csi"),emit:vcf
    tuple val(meta),path("*.md5")
    path("versions.yml"),emit:versions
script:
    def prefix=task.ext.prefix?:"${meta.id}.cnvnator"
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP
find VCFS/ \\( -name "*.bcf" -o -name "*.vcf.gz" \\) > TMP/jeter.list
bcftools concat -a -D --threads ${task.cpus} --write-index --file-list TMP/jeter.list -O b9 -o TMP/jeter.bcf

mv TMP/jeter.bcf ${prefix}.bcf
mv TMP/jeter.bcf.csi ${prefix}.bcf.csi

md5sum ${prefix}.bcf  > ${prefix}.bcf.md5

cat << END_VERSIONS > versions.yml
${task.process}:
    bcftools: \$(bcftools version | awk '(NR==1)  {print \$NF}')
END_VERSIONS
"""

stub:
 def prefix=task.ext.prefix?:meta.id
"""
touch versions.yml ${prefix}.bcf ${prefix}.bcf.csi  ${prefix}.bcf.md5
"""
}

process CNVNATOR_CONCAT_BEDS {
tag "${meta.id}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
    tuple val(meta),path("BEDS/*")
output:
    tuple val(meta),path("*.bed.gz"),path("*.bed.gz.tbi"),emit:vcf
    path("versions.yml"),emit:versions
script:
    def prefix=task.ext.prefix?:meta.id
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP
export LC_ALL=C

echo "#chrom start end sample CNV_type coordinates CNV_size normalized_RD e-val1 e-val2 e-val3 e-val4 q0" | tr " " "\t" >  TMP/jeter.bed

gunzip -c  BEDS/*.bed.gz | sort  -S ${task.memory.kilo} -T TMP -t '\t'  -k1,1 -k2,2n >> TMP/jeter.bed
bgzip TMP/jeter.bed
tabix -f -p bed TMP/jeter.bed.gz

mv TMP/jeter.bed.gz ${prefix}.bed.gz
mv TMP/jeter.bed.gz.tbi ${prefix}.bed.gz.tbi

cat << END_VERSIONS > versions.yml
${task.process}:
    tabix: \$(tabix --version | awk '(NR==1)  {print \$NF}')
END_VERSIONS
"""

stub:
 def prefix=task.ext.prefix?:meta.id
"""
touch versions.yml ${prefix}.bed.gz ${prefix}.bed.gz.tbi
"""
}
