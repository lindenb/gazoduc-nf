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
include { BWA_INDEX       } from '../../modules/bwa/index/main.nf'
include { SPADES_ASSEMBLY } from '../../modules/spades/assembly'
include { parseBoolean    } from '../../modules/utils/functions.nf'

workflow UNMAPPED {
take:
    metadata
    acn_file
    bams // tuple (meta,bam,bai,fasta,fai)
main:
    versions  = Channel.empty()
    multiqc  = Channel.empty()
    FETCH_ACNS(acn_file)
    versions  = versions.mix(FETCH_ACNS.out.versions)

    


    BWA_INDEX(FETCH_ACNS.out.fasta)
    //versions  = versions.mix(BWA_INDEX.out.versions)

    MAP_UNMAPPED(
        BWA_INDEX.out.bwa_index,
        bams
        )
    versions  = versions.mix(MAP_UNMAPPED.out.versions)

    MERGE_CONTAMINANTS_MAPPING(
        FETCH_ACNS.out.database,
        MAP_UNMAPPED.out.tsv
            .map{meta,tsv->tsv}
            .collect()
            .map{[[id:"contaminations"],it.sort()]}
        )
    versions  = versions.mix(MERGE_CONTAMINANTS_MAPPING.out.versions)

    SPADES_ASSEMBLY(
        MAP_UNMAPPED.out.fastq.map{[it[0],it[1],[]/* no fastq_2 */]}
        )
    versions  = versions.mix(SPADES_ASSEMBLY.out.versions)


    
    MERGE_SPADES(
        SPADES_ASSEMBLY.out.assembled_contigs_tsv
            .map{meta,tsv->tsv}
            .collect()
            .map{[[id:"contaminations"],it.sort()]}
        )
    versions  = versions.mix(MERGE_SPADES.out.versions)
emit:
    versions
    multiqc
    spades_tsv_gz = MERGE_SPADES.out.tsv_gz

}

process FETCH_ACNS {
    tag "${acns_file.name}"
    label "label_single"
    conda "${moduleDir}/../../conda/bioinfo.01.yml"
    input:
        tuple val(meta),path(acns_file)
    output:
        tuple val(meta),path("*.fa"),emit:fasta
        tuple val(meta),path("*.fa.fai"),emit:fai
        tuple val(meta),path("*.dict"),emit:dict
        tuple val(meta),path("*.tsv"),emit:database
        path("versions.yml"),emit:versions
    script:
        def prefix = task.ext.prefix?:acns_file.baseName+".db"
"""
mkdir -p TMP
mkdir -p BWAINDEX

curl --version 1>&2
env | grep -i proxy

set -x
cat "${acns_file}"|\\
    grep -v '^#' |\\
    tr -s " " | tr " " "\n" | sort | uniq |\\
    grep -v '^\$'  > TMP/jeter.txt

test -s TMP/jeter.txt

cat TMP/jeter.txt |\\
    xargs -n 20  echo  | tr " " "," |\\
    while read ACN
    do
        if test -f TMP/jeter.fa
        then
            # prevent too many calls to NCBI
            sleep 10
        fi

        curl -k -L -o TMP/jeter0.fa "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=\${ACN}&rettype=fasta" 
        grep -v '^\$' TMP/jeter0.fa >> TMP/jeter.fa
        rm TMP/jeter0.fa

    done

grep '^>' TMP/jeter.fa |\\
    cut -c 2- |\\
    awk '{printf("%s\t%s\\n",\$1,\$0);}' |\\
    sort -t '\t' -T TMP -k1,1 > ${prefix}.tsv

cut -d ' ' -f 1  TMP/jeter.fa > ${prefix}.fa
samtools faidx ${prefix}.fa
samtools dict ${prefix}.fa > ${prefix}.dict


cat << END_VERSIONS > versions.yml
${task.process}:
    samtools : \$(samtools version | awk 'NR==1 {print \$NF}')
END_VERSIONS
"""
stub:
        def prefix = task.ext.prefix?:acns_file.baseName+".db"

"""
touch ${prefix}.fa ${prefix}.fa.fai ${prefix}.dict ${prefix}.tsv
touch versions.yml
"""
}



process MAP_UNMAPPED {
    tag "${meta.id}"
    label "label_short"
    conda "${moduleDir}/../../conda/bwa.yml"
    afterScript "rm -rf TMP"
    input:
        tuple val(meta5),path(bwaDir)
        tuple val(meta),path(bam),path(bai),path(fasta),path(fai),path(dict)
    output:
        tuple val(meta), path("*.contaminations.tsv"),emit:tsv
        tuple val(meta), path("*.fastq.gz"),optional:true,emit:fastq // still unmapped, for spades
        path("versions.yml"),emit:versions
    script:
        def prefix=task.ext.prefix?:meta.id
        def cpusargs1 = task.cpus > 5? "--threads ${task.cpus -5 }":""
        def tresholdN = 0.3
        def repeats=6
	    def min_size  = task.ext.min_size?:50
        def fast = parseBoolean(task.ext.fast?:false) /*      Output the unmapped reads at the end of the file */
"""
mkdir -p TMP
set -x

cat << 'EOF' > TMP/jeter.awk
    {
    OFS="\t";
    S=\$2;
    L1=length(S)*1.0;
    if(L1< ${min_size} ) next;
    
    S=\$2;
    gsub(/[ATat]/,"",S);
    L2=length(S)*1.0;
    if(L2/L1 <= 0.3 || L2/L1>=0.7 ) next; 

    S=\$2; 
    gsub(/[ATGCatgc]/,"",S);
    L2=length(S)*1.0;
    if(L2/L1 >= ${tresholdN} ) next;  
    
    S=\$2; 
    gsub(/A{${repeats},}/,"",S);
    L2=length(S)*1.0;
    if(L2/L1 <= ${tresholdN} ) next; 
    
    S=\$2; 
    gsub(/C{${repeats},}/,"",S);
    L2=length(S)*1.0;
    if(L2/L1 <= ${tresholdN} ) next;
    
    S=\$2; 
    gsub(/G{${repeats},}/,"",S);
    L2=length(S)*1.0;
    if(L2/L1 <= ${tresholdN} ) next;

    S=\$2; 
    gsub(/T{${repeats},}/,"",S);
    L2=length(S)*1.0;
    if(L2/L1 <= ${tresholdN} ) next; 
    
    S=\$2; 
    gsub(/CA{${repeats},}/,"",S);
    L2=length(S)*1.0;
    if(L2/L1 <= ${tresholdN} ) next; 

    S=\$2; 
    gsub(/TG{${repeats},}/,"",S);
    L2=length(S)*1.0;
    if(L2/L1 <= ${tresholdN} ) next; 

    printf("%s\\n%s\\n%s\\n%s\\n",\$1,\$2,\$3,\$4);
    }

EOF

samtools view ${cpusargs1} --uncompressed --reference ${fasta} -f 4 -F 3840 "${bam}" ${fast?"\"*\"":""} |\\
    samtools fastq -n - |\\
    paste - - - - |\\
    awk -F '\t' -f TMP/jeter.awk |\\
    gzip > TMP/unmapped.fastq.gz


bwa mem \\
     -R '@RG\\tID:${meta.id}\\tSM:${meta.id}\\tLB:${meta.id}\\tCN:NantesBird\\tPL:ILLUMINA' \\
    -t ${task.cpus} \\
    `find ${bwaDir}/ -name "*.amb" | sed 's/\\.amb\$//'` \\
    TMP/unmapped.fastq.gz |\\
    samtools view  -o TMP/jeter.bam -

# still unmapped
samtools view --uncompressed --reference ${fasta} -f 4 -F 3840 TMP/jeter.bam |\\
    samtools fastq -n - |\\
    gzip  > TMP/unmapped2.fastq.gz

samtools view -F 3844  --uncompressed TMP/jeter.bam |\\
samtools sort  -m '${task.memory.giga}G' --threads '${task.cpus}' -o TMP/jeter2.bam -O BAM -T TMP/tmp -
samtools index -@ ${task.cpus} "TMP/jeter2.bam"

samtools idxstat TMP/jeter2.bam |\\
    awk -F '\t' '(\$3!=0 && \$1!="*") {printf("%s\t%s\t${meta.id}\\n",\$1,\$3);}' |\
    sort -T TMP -k1,1 > TMP/jeter.tsv

mv TMP/jeter.tsv ${prefix}.contaminations.tsv

if [[ \$(gunzip -c TMP/unmapped2.fastq.gz | wc -l ) -gt 3 ]]
then
    mv TMP/unmapped2.fastq.gz ${prefix}.unmapped.fastq.gz
fi

cat << END_VERSIONS > versions.yml
${task.process}:
    samtools : \$(samtools version | awk 'NR==1 {print \$NF}')
    bwa : \$(bwa 2>&1 | grep Version | tr -d ":")
END_VERSIONS
"""

stub:
        def prefix=task.ext.prefix?:meta.id
"""
touch versions.yml ${prefix}.contaminations.tsv  ${prefix}.fastq.gz
"""
}


process MERGE_CONTAMINANTS_MAPPING {
tag "${meta.id}"
label "label_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
    tuple val(meta1),path("database.tsv")
    tuple val(meta ),path("TSV/*")
output:
    tuple val(meta),path("*.tsv",arity:"1"),emit:tsv
    path("versions.yml"),emit:versions
script:
    def prefix=task.ext.prefix?:"contaminations"
"""
mkdir -p TMP

find TSV/ -type l -exec cat '{}' ';' |\\
    sort -t '\t' -k1,1 -T TMP > TMP/jeter.a

join -t '\t' -1 1 -2 1 -o '1.1,1.2,1.3,2.2' TMP/jeter.a database.tsv |\
    sort -t '\t' -T TMP -k2,2n > ${prefix}.tsv

cat << END_VERSIONS > versions.yml
${task.process}:
    join : todo
    sort : todo
END_VERSIONS
"""

stub:
   def prefix=task.ext.prefix?:"contaminations"
"""
touch versions.yml ${prefix}.tsv
"""
}

process MERGE_SPADES {
tag "${meta.id}"
label "label_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
    tuple val(meta ),path("TSV/*")
output:
    tuple val(meta),path("*.tsv.gz",arity:"1"),emit:tsv_gz
    path("versions.yml"),emit:versions
script:
    def prefix=task.ext.prefix?:"spades"
"""
mkdir -p TMP

echo -e 'LENGTH\tNAME\tSAMPLE\tFASTA' >> TMP/jeter.a

find TSV/ -type l -exec cat '{}' ';' |\\
    gunzip -c |\\
    sort -t '\t' -k1,1nr -T TMP  >> TMP/jeter.a

gzip --best TMP/jeter.a
mv  TMP/jeter.a.gz ${prefix}.tsv.gz


cat << END_VERSIONS > versions.yml
${task.process}:
    sort : todo
END_VERSIONS
"""
stub:
   def prefix=task.ext.prefix?:"spades"
"""
touch versions.yml ${prefix}.tsv.gz
"""
}
