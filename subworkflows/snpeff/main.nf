include {k1_signature} from '../../modules/utils/k1.nf'

workflow SNPEFF {
    take:
        meta
        fasta
        fai
        dict
        vcf
    main:
        DOWNLOAD(fai)
        ANNOTATE(DOWNLOAD.out.db, DOWNLOAD.out.database_name, vcf)
    emit:
        vcf = ANNOTATE.out.vcf
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
script:
    def k1= k1_signature()
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

snpEff -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  download -dataDir  "\${PWD}/SNPEFFX"  "\${DB}"

test -s SNPEFFX/*/snpEffectPredictor.bin
mv SNPEFFX "SNPEFF"

mv TMP/jeter.database "\${DB}.database"
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
	tuple val(meta),path("*.bcf"),path("*.bcf.csi"),emit:vcf
script:
    def prefix=task.ext.prefix?:vcf.simpleName+".snpeff"
"""
hostname 1>&2
mkdir -p TMP OUTPUT

set -o pipefail

bcftools view '${vcf}' -O v |\\
	snpEff -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP eff \\
        -dataDir "\${PWD}/${config}" \\
		-nodownload \\
        -noLog \\
        -noStats \\
        -lof \\
        `cat ${db_name}` > TMP/jeter1.vcf

bcftools sort --max-mem '${task.memory.giga}G' -T TMP/tmp -O b -o TMP/${prefix}.bcf TMP/jeter1.vcf
bcftools index --threads ${task.cpus} TMP/${prefix}.bcf

mv TMP/${prefix}.bcf ./
mv TMP/${prefix}.bcf.csi ./
"""
}
