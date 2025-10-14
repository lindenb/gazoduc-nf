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

process HAPLOTYPECALLER {
tag "${meta.id?:""} ${optional_bed?optional_bed.name:""}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.02.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta4),path("REFS/*") //bag fasta/fai will be used change the dict if ref(bam) is not the fasta
    tuple val(meta ),path(bam),path(bai),path(optional_bed)
output:
    tuple val(meta),path("*.g.vcf.gz"),path("*.g.vcf.gz.tbi"),path(optional_bed),emit:gvcf
    path("versions.yml"),emit:versions
script:
   def prefix0 = (meta.id?:"${bam.name}")+(optional_bed?"."+optional_bed.baseName:"")
   def prefix= task.ext.prefix?:prefix0
   def args1 = task.ext.args1?:"-G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation"
   def args2 =task.ext.args2?:""
   def jvm = task.ext.jvm?:"-XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -DGATK_STACKTRACE_ON_USER_EXCEPTION=true"
"""
hostname 1>&2
mkdir -p TMP
# if not build because no other ref
mkdir -p REFS
find \${PWD}/REFS/ \\( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \\) >> TMP/references.txt
find REFS/ 1>&2

samtools samples -f "${fasta}" -F TMP/references.txt "${bam}" | cut -f1,3 | head -n1 | while read SAMPLE REF
    do

        # test fasta is known
        test "\${REF}" != "."

        
        if ${optional_bed?true:false}
        then
            if cmp "\${REF}.fai" "${fai}"
            then
               cp -v "${optional_bed}" TMP/${optional_bed.name}
            else
                # no the same reference ? change the BED according to chr notation
                jvarkit ${jvm} bedrenamechr \\
                        -f "\${REF}" --column 1 --convert SKIP "${optional_bed}" > TMP/${optional_bed.name}
                
                # never empty file
                if ! test -s TMP/${optional_bed.name}
                then
                        awk -F '\t' '{printf("%s\t0\t1\\n",\$1);}' "\${REF}.fai" | tail -n 1 > TMP/${optional_bed.name}
                fi

            fi
        fi

        gatk --java-options "${jvm}" HaplotypeCaller \\
            ${optional_bed?"-L TMP/${optional_bed.name}":""} \\
            -R "\${REF}" \\
            -I "${bam}" \\
            -ERC GVCF \\
            ${args1} ${args2}\\
            -O "TMP/jeter.g.vcf.gz"


        # rename sample if needed
        if test "${meta.id}" != "\${SAMPLE}"
        then
                gatk --java-options "${jvm}"  RenameSampleInVcf \
                        -INPUT "TMP/jeter.g.vcf.gz" \
                        -OUTPUT TMP/jeter2.g.vcf.gz \
                        -NEW_SAMPLE_NAME "${meta.id}"

                rm -f TMP/jeter.g.vcf.gz.tbi
                mv TMP/jeter2.g.vcf.gz "TMP/jeter.g.vcf.gz"
        fi


        # not the original fasta ?
        if ! cmp "\${REF}.fai" "${fai}"
        then
            jvarkit ${jvm} vcfsetdict \\
                -n SKIP \\
                -R "${fasta}" \\
                "TMP/jeter.g.vcf.gz" > TMP/jeter2.vcf

            rm -f TMP/jeter.g.vcf.gz.tbi
            bcftools sort  -T TMP/sort  --max-mem "${task.memory.giga}G" \\
                -O z -o "TMP/jeter.g.vcf.gz" TMP/jeter2.vcf
        fi

    done

bcftools index --threads ${task.cpus} --force --tbi "TMP/jeter.g.vcf.gz"

mv "TMP/jeter.g.vcf.gz" ${prefix}.g.vcf.gz
mv "TMP/jeter.g.vcf.gz.tbi" ${prefix}.g.vcf.gz.tbi

cat << EOF > versions.yml
${task.process}:
    gatk: "\$( gatk --version 2> /dev/null  | paste -s -d ' ' )"
EOF
"""


stub:
   def prefix = (meta.id?:"${bam.name}")+(optional_bed?"."+optional_bed.baseName:"")

"""
touch versions.yml ${prefix}.g.vcf.gz ${prefix}.g.vcf.gz.tbi
"""
}
