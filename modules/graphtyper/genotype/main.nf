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
process GRAPHTYPER {
tag "${bed.name}"
label "process_single"
afterScript "rm -rf TMP TMP2"
conda "${moduleDir}/../../../conda/graphtyper.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta4),path(optional_covlen)
    tuple val(meta),path("BAMS/*"),path(bed)
output:
    tuple val(meta),path("*.bcf",arity:'1'),path("*.csi",arity:'1'),emit:vcf
    path("versions.yml"),emit:versions
script:
    def args = task.ext.args?:""
    def prefix = task.ext.prefix?:"\${MD5}"
"""

mkdir -p TMP
mkdir -p TMP2
find BAMS/ -name "*am" | sort -V > TMP/bams.list


if ${optional_covlen?true:false}
then
    # join samples name and content of optional_covlen
    samtools samples < TMP/bams.list | sort -t '\t' -T TMP -k1,1 --unique | cut -f1,2 > TMP/jeter1.tsv
    sort -t '\t' -T TMP -k1,1 --unique "${optional_covlen}" | cut -f1,2 > TMP/jeter2.tsv

    # paranoid
    cut -f 1 TMP/jeter1.tsv > TMP/jeter.a
    cut -f 1 TMP/jeter2.tsv > TMP/jeter.b
    # bam without data ?
    comm -23 TMP/jeter.a TMP/jeter.b > TMP/jeter.dup && test ! -s  TMP/jeter.dup

    join -t '\t' -1 1 -2 1 -o 1.2 TMP/jeter1.tsv TMP/jeter2.tsv > TMP/bams.list
    join -t '\t' -1 1 -2 1 -o 2.2 TMP/jeter1.tsv TMP/jeter2.tsv > TMP/cov.length.tsv

    test -s TMP/cov.length.tsv
fi

test -s  TMP/bams.list

MD5=`echo ${bed}  | sha1sum | cut -d ' ' -f1`

awk -F '\t' '{printf("%s:%d-%s\\n",\$1,int(\$2)+1,\$3);}' "${bed}"  > TMP/jeter.intervals

graphtyper genotype \\
        "${fasta}" \\
        ${optional_covlen?"--avg_cov_by_readlen=TMP/cov.length.tsv":""} \\
        --output=TMP2 \\
        --force_no_copy_reference \\
        --force_use_input_ref_for_cram_reading \\
        --sams=TMP/bams.list \\
        --region_file=TMP/jeter.intervals \\
        --threads=${task.cpus} \\
        ${args}

rm -rf "TMP2/input_sites"


rm -f TMP/jeter.list
find TMP2 -type f -name "*.vcf.gz"  | while read F
do
        bcftools query -l "\${F}" | sort > TMP/ordered.samples.txt
        bcftools view  --threads ${task.cpus} -O b -o "\${F}.bcf"  --samples-file TMP/ordered.samples.txt "\${F}"
        bcftools index --threads ${task.cpus} "\${F}.bcf"
	    
        echo "\${F}.bcf" >> TMP/jeter.list
        rm TMP/ordered.samples.txt
        rm "\${F}" "\${F}.tbi"
done


test -s TMP/jeter.list

# merge all
bcftools concat --threads ${task.cpus} --allow-overlaps --remove-duplicates -Ob -o "${prefix}.bcf" --file-list TMP/jeter.list
bcftools index --threads ${task.cpus} ${prefix}.bcf


cat << END_VERSIONS > versions.yml
${task.process}:
    graphtyper: \$(graphtyper --help | tail -n 1 | sed 's/^   //')
END_VERSIONS
"""
}
