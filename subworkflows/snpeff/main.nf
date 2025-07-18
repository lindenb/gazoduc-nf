include {k1_signature} from '../../modules/utils/k1.nf'

workflow SNPEFF {
    take:
        meta
        fasta
        fai
        dict
        vcf
    main:
        versions = Channel.empty()
        DOWNLOAD(fai)
        versions = versions.mix(DOWNLOAD.out.versions)

        ANNOTATE(DOWNLOAD.out.db, DOWNLOAD.out.database_name, vcf)
        versions = versions.mix(ANNOTATE.out.versions)
    emit:
        database = DOWNLOAD.out.db
        database_name = DOWNLOAD.out.database_name
        vcf = ANNOTATE.out.vcf
        versions
}


process DOWNLOAD {
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
    tuple val(meta),path(fai)
output:
    tuple val(meta),path("SNPEFF"),emit:db
    tuple val(meta),path("*.database"),emit:database_name
    path("versions.yml"),emit:versions
script:
    def k1= k1_signature()
    def snpeff_database_directory = task.ext.snpeff_database_directory?:"NODIR"
"""
mkdir -p SNPEFFX TMP

cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg38}\tGRCh38.99
1:${k1.hg19}\tGRCh37.75
EOF


awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sed 's/^chr//' | sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv
join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv | sort | uniq > TMP/jeter.database

DB=`cat TMP/jeter.database`

test ! -z "\${DB}"

if test -f "${snpeff_database_directory}/\${DB}/snpEffectPredictor.bin"
then

    ln -s "${snpeff_database_directory}" "SNPEFF"

else

snpEff  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  download -dataDir  "\${PWD}/SNPEFFX"  "\${DB}"

test -s SNPEFFX/*/snpEffectPredictor.bin
mv SNPEFFX "SNPEFF"
fi

mv TMP/jeter.database "\${DB}.database"


cat << END_VERSIONS > versions.yml
"${task.process}":
	snpeff: "todo"
END_VERSIONS
"""
}


process ANNOTATE {
tag "${meta.id}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"

input:
	tuple val(meta1),path(config)
    tuple val(meta2),path(db_name)
	tuple val(meta ),path(vcf),path(vcfidx)
output:
	tuple val(meta),path("*.bcf"),path("*.bcf.csi"),optional:true,emit:vcf
    path("versions.yml"),emit:versions
script:
    def prefix=task.ext.prefix?:vcf.baseName+".snpeff"
    def accessions = task.ext.accessions?:"" // keep blank, no filter by default. Otherwise you want SO:0001818,SO:0001629"
    def args1 = task.ext.args1?:"-noLog -noStats -lof"
    def args2 = task.ext.args2?:""

"""
hostname 1>&2
mkdir -p TMP OUTPUT

set -o pipefail

bcftools view '${vcf}' -O v |\\
	snpEff -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP eff \\
        -dataDir "\${PWD}/${config}" \\
		-nodownload \\
        ${args1} \\
        `cat ${db_name}` > TMP/jeter1.vcf

if ${!accessions.trim().isEmpty()}
then
    jvarkit -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP \\
        vcffilterso \\
        ${args2} \\
        --acn "${accessions}" \\
        TMP/jeter1.vcf >  TMP/jeter2.vcf
    
    mv TMP/jeter2.vcf TMP/jeter1.vcf
fi

bcftools sort --max-mem '${task.memory.giga}G' -T TMP/tmp -O b -o TMP/jeter.bcf TMP/jeter1.vcf
bcftools index --threads ${task.cpus} TMP/jeter.bcf


if test \$(bcftools index -s TMP/jeter.bcf |wc -l) -gt 0
then
    mv TMP/jeter.bcf ${prefix}.bcf 
    mv TMP/jeter.bcf.csi ${prefix}.bcf.csi
fi


cat << END_VERSIONS > versions.yml
"${task.process}":
	snpeff: "todo"
END_VERSIONS
"""
}
