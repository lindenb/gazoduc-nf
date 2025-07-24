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
process SPADES_ASSEMBLY {
    tag "${meta.id}"
    label "process_medium"
    afterScript "rm -rf TMP"
    conda "${moduleDir}/../../../conda/spades.yml"
    input:
        tuple val(meta),path(fastq_1),path(optional_fastq_2)
    output:
        tuple val(meta),path("*.assembled_contigs.fa.gz" ),optional:true,emit:assembled_contigs_fa
        tuple val(meta),path("*.assembled_contigs.tsv.gz"),optional:true,emit:assembled_contigs_tsv
        path("versions.yml"),emit:versions
    script:
        def has_fq2 = optional_fastq_2?true:false
        def args1 = task.ext.args1?:" -k 21,33,55,77,99,127 --careful --sc "
        def prefix = task.ext.prefix?:meta.id+".spades"
    """
    mkdir -p TMP
    export TMPDIR=\${PWD}/TMP

    spades.py ${args1} \\
        -t ${task.memory.giga} \\
        -m ${task.cpus} \\
        -o TMP/jeter_assembly \\
        ${has_fq2 ? "-1 \"${fastq_1}\" -2 \"${optional_fastq_2}\" ":"-s \"${fastq_1}\""}


    if test -s TMP/jeter_assembly/misc/assembled_contigs.fasta
    then
        

        awk '/^>/ {printf("%s%s\t%s\t",(N>0?"\\n":""),substr(\$0,2),"${meta.id}");N++;next;} {printf("%s",\$0);} END {printf("\\n");}' TMP/jeter_assembly/misc/assembled_contigs.fasta |\\
        awk -F '\t' '{printf("%d\t%s",length(\$3),\$0);}' |\\
        gzip --best >  TMP/jeter_assembly/misc/assembled_contigs.tsv

        gzip --best TMP/jeter_assembly/misc/assembled_contigs.fasta

        mv "TMP/jeter_assembly/misc/assembled_contigs.fasta.gz" "${prefix}.assembled_contigs.fa.gz"
        mv "TMP/jeter_assembly/misc/assembled_contigs.tsv.gz" "${prefix}.assembled_contigs.tsv.gz"
    fi

cat << EOF > versions.yml
${task.process}:
	spades: \$(spades.py --version 2>&1 | sed -n 's/^.*SPAdes genome assembler v//p')
EOF
    """

}