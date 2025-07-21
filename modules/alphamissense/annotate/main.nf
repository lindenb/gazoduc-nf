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

process ALPHAMISSENSE_ANNOTATE {
tag "${vcf.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(tabix),path(tbi),path(header)
    tuple val(meta ),path(vcf),path(vcf_idx)
output:
    tuple val(meta),path("*.bcf"),path("*.csi"),emit:vcf
	path("versions.yml"),emit:versions
script:
   	def TAG = task.ext.tag?:"ALPHAMISSENSE"
	def prefix = task.ext.prefix?:vcf.baseName+".alphamissense"
"""
hostname 1>&2
mkdir -p TMP OUTPUT

bcftools annotate \\
    --threads ${task.cpus} \\
    -a "${tabix}" \\
    -h "${header}" \\
    -c "CHROM,POS,REF,ALT,${TAG}_PATHOGENOCITY,${TAG}_CLASS" \\
    --merge-logic "${TAG}_PATHOGENOCITY:max,${TAG}_CLASS:unique" \\
    -O b -o TMP/${prefix}.bcf '${vcf}'
 
bcftools index \\
    --threads ${task.cpus} \\
    -f TMP/${prefix}.bcf

mv TMP/*.bcf ./
mv TMP/*.bcf.csi ./

cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
}