include { BWA_INDEX       } from '../../modules/bwa/index/main.nf'
include { SPADES_ASSEMBLY } from '../../modules/spades/assembly'
workflow CONTAMINATIONS {
take:
    meta
    acn_file
    bams // tuple (meta,bam,fasta,fai)
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

    SPADES_ASSEMBLY(
        MAP_UNMAPPED.out.fastq.map{[it[0],it[1],[]/* no fastq_2 */]}
        )
    ersions  = versions.mix(SPADES_ASSEMBLY.out.versions)

    MERGE(
        FETCH_ACNS.out.database,
        MAP_UNMAPPED.out.tsv
            .map{it[1]}
            .collect()
            .map{[[id:"contaminations"],it]}
        )
    
    MERGE_SPADES(
        SPADES_ASSEMBLY.out.assembled_contigs_tsv
            .map{it[1]}
            .collect()
            .map{[[id:"contaminations"],it]}
        )
    versions  = versions.mix(MERGE_SPADES.out.versions)
emit:
    versions
    multiqc
    output = MERGE.out.output

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
}



process MAP_UNMAPPED {
    tag "${meta.id}"
    label "label_single"
    conda "${moduleDir}/../../conda/bwa.yml"
    afterScript "rm -rf TMP"
    input:
        tuple val(meta5),path(bwaDir)
        tuple val(meta),path(bam),path(fasta),path(fai)
    output:
        tuple val(meta), path("*.contaminations.tsv"),emit:tsv
        tuple val(meta), path("*.fastq.gz"),emit:fastq // still unmapped, for spades
        path("versions.yml"),emit:versions
    script:
        def prefix=task.ext.prefix?:meta.id
        def tresholdN = 0.5
        def repeats=6
"""
mkdir -p TMP
set -x

cat << 'EOF' > TMP/jeter.awk
    {
    OFS="\t";
    S=\$2;
    L1=length(S)*1.0;
    if(L1<50) next;
    
    S=\$2;
    gsub(/[ATat]/,"",S);
    L2=length(S)*1.0;
    if(L2/L1 <= 0.2 || L2/L1>=0.8 ) next; 

    S=\$2; 
    gsub(/[^ATGCatgc]/,"",S);
    L2=length(S)*1.0;
    if(L2/L1 <= ${tresholdN} ) next;  
    
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

samtools view --uncompressed --reference ${fasta} -f 4 -F 3840 "${bam}" |\\
    samtools fastq -n - |\\
    paste - - - - |\\
    awk -F '\t' -f TMP/jeter.awk |\\
    gzip --best > TMP/unmapped.fastq.gz


bwa mem \\
     -R '@RG\\tID:${meta.id}\\tSM:${meta.id}\\tLB:${meta.id}\\tCN:NantesBird\\tPL:ILLUMINA' \\
    -t ${task.cpus} \\
    `find ${bwaDir}/ -name "*.amb" | sed 's/\\.amb\$//'` \\
    TMP/unmapped.fastq.gz |\\
    samtools view --uncompressed -o TMP/jeter.bam -

# still unmapped
samtools view --uncompressed --reference ${fasta} -f 4 -F 3840 TMP/jeter.bam |\\
    samtools fastq -n - |\\
    gzip --best > TMP/unmapped2.fastq.gz

samtools view -F 3844  --uncompressed TMP/jeter.bam |\\
samtools sort  -m '${task.memory.giga}G' --threads '${task.cpus}' -o TMP/jeter2.bam -O BAM -T TMP/tmp -
samtools index -@ ${task.cpus} "TMP/jeter2.bam"

samtools idxstat TMP/jeter2.bam |\\
    awk -F '\t' '(\$3!=0 && \$1!="*") {printf("%s\t%s\t${meta.id}\\n",\$1,\$3);}' |\
    sort -T TMP -k1,1 > TMP/jeter.tsv

mv TMP/jeter.tsv ${prefix}.contaminations.tsv

mv TMP/unmapped2.fastq.gz ${prefix}.unmapped.fastq.gz

cat << END_VERSIONS > versions.yml
${task.process}:
    samtools : \$(samtools version | awk 'NR==1 {print \$NF}')
    bwa : \$(bwa 2>&1 | grep Version | tr -d ":")
END_VERSIONS
"""
}


process MERGE {
tag "${meta.id}"
label "label_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
    tuple val(meta1),path("database.tsv")
    tuple val(meta ),path("TSV/*")
output:
    tuple val(meta),path("*.tsv",arity:"1"),emit:output
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
}

process MERGE_SPADES {
tag "${meta.id}"
label "label_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
    tuple val(meta ),path("TSV/*")
output:
    tuple val(meta),path("*.tsv",arity:"1"),emit:output
    path("versions.yml"),emit:versions
script:
    def prefix=task.ext.prefix?:"spades"
"""
mkdir -p TMP

echo -e 'LENGTH\tNAME\tSAMPLE\tFASTA' >> TMP/jeter.a

find TSV/ -type l -exec cat '{}' ';' |\\
    sort -t '\t' -k1,1nr -T TMP >> TMP/jeter.a


mv TMP/jeter.a ${prefix}.tsv


cat << END_VERSIONS > versions.yml
${task.process}:
    sort : todo
END_VERSIONS
"""
}
