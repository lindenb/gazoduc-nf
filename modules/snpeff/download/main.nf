
process SNPEFF_DOWNLOAD {
label "process_single"
tag "${snpeff_db}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	val(snpeff_db)
output:
        path("SNPEFF.${snpeff_db}"),emit:output
        paht("versions.yml"),emit:versions
script:
"""
mkdir -p SNPEFFX TMP
snpEff -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  download -dataDir  "\${PWD}/SNPEFFX"  '${snpeff_db}'

test -s SNPEFFX/*/snpEffectPredictor.bin
mv SNPEFFX "SNPEFF.${snpeff_db}"

cat << END_VERSIONS > versions.yml
"${task.process}":
    db: todo
END_VERSIONS
"""
}



