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
process GATHER_VCFS {
tag "${meta.id?:""}"
label "process_short"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	tuple val(meta),path("VCFS/*")
output:
	tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
    path("versions.yml"),emit:versions
when:
    task.ext.when == null || task.ext.when
script:
    def prefix = task.ext.prefix?:"\${MD5}.gathervcfs"
"""
hostname 1>&2
mkdir -p TMP
# gonna fail if there is a BCF, paranoid
find VCFS/ \\( -name "*.vcf.gz" -o -name "*.bcf" \\) > TMP/jeter.list
MD5=`cat TMP/jeter.list | md5sum | cut -d ' ' -f1`

gatk  --java-options "-Xmx${task.memory.giga}g -XX:-UsePerfData -Djava.io.tmpdir=TMP" GatherVcfs  \\
	--INPUT TMP/jeter.list \\
	--REORDER_INPUT_BY_FIRST_VARIANT \\
	--OUTPUT TMP/jeter.vcf.gz	

bcftools index --force --threads "${task.cpus}" TMP/jeter.vcf.gz

mv TMP/jeter.vcf.gz "${prefix}.vcf.gz"
mv TMP/jeter.vcf.gz.tbi "${prefix}.vcf.gz.tbi"

cat << END_VERSIONS > versions.yml
"${task.process}":
	gatk: todo
	bcftools: todo
END_VERSIONS
"""

stub:
"""
touch versions.yml "${meta.id}.gb.vcf.gz" "${meta.id}.gb.vcf.gz.tbi"
"""
}