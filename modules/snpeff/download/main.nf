
process SNPEFF_DOWNLOAD {
label "process_single"
tag "${snpeff_db}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP SNPEFFX"
input:
    tuple val(meta),path(fai)
output:
    tuple val(meta),path("*.DB"),path("*.name.txt"),emit:database
    path("versions.yml"),emit:versions
script:
    def snpeff_database_directory = task.ext.snpeff_database_directory?:"NODIR"
    def snpeff_db = task.ext.db?:""
    if(snpeff_db.isEmpty() && meta.ucsc_name && meta.ucsc_name.equals("hg19")) {
        snpeff_db = "GRCh37.75"
        }
    else if(snpeff_db.isEmpty() && meta.ucsc_name && meta.ucsc_name.equals("hg38")) {
        snpeff_db = "GRCh38.99"
        }
    else if(snpeff_db.isEmpty())
        {
        throw new IllegalArgumentException("${task.process} undefined snpeff_db");
        }
    def prefix = task.ext.prefix?:"SNPEFF"
"""
mkdir -p SNPEFFX TMP

test ! -z "${snpeff_db}"

if test -f "${snpeff_database_directory}/${snpeff_db}/snpEffectPredictor.bin"
then

    ln -s "${snpeff_database_directory}" "${prefix}.${snpeff_db}.DB"

else

    snpEff \${JAVA_PROXY:-} -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  \\
        download -dataDir  "\${PWD}/SNPEFFX"  "${snpeff_db}"

    test -s SNPEFFX/*/snpEffectPredictor.bin
    mv SNPEFFX "${prefix}.${snpeff_db}.DB"
fi

echo '${snpeff_db}' > "${snpeff_db}.name.txt"


cat << END_VERSIONS > versions.yml
"${task.process}":
        snpeff: "todo"
        snpeffdb: "${snpeff_db}"
END_VERSIONS
"""
}



