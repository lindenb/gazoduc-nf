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
process SNPEFF_DOWNLOAD {
label "process_single"
tag "${meta.id}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP SNPEFFX"
input:
    tuple val(meta),path(snpeff_database_directory) //load directoy, may not exit but cannot be null
output:
    tuple val(meta),path("*.DB"),path("*.name.txt"),emit:database
    path("versions.yml"),emit:versions
script:
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

stub:
   def snpeff_db = "GRCh37.75"
   def prefix = task.ext.prefix?:"SNPEFF"
"""
mkdir -p  "${prefix}.${snpeff_db}.DB"
touch versions.yml "${snpeff_db}.name.txt"
"""
}



