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
process REGENIE_STEP1 {
label "process_single"
label "process_high_memory"
tag "${meta1.id}"
conda "${moduleDir}/../../conda/regenie.yml"
afterScript "rm -rf TMP"
input:
        tuple val(meta1),path(bgen_files)
        tuple val(meta2),path(covariates)
        tuple val(meta3),path(pheno_files)
        tuple val(meta4),path(plink_files)
output:
        tuple val(meta1),path("*_pred.list"),emit:output
        tuple val(meta1),path("*.log"),emit:log
        path("versions.yml"),emit:versions
script:
        def prefix = task.ext.prefix?:"{meta1.id}step1"
        def pgen = bgen_files.find{it.name.endsWith(".pgen")}
        def keep_rs = plink_files.find{it.name.endsWith("keep.id.txt")}
        def args = task.ext.args?:"--bsize 1000 --bt --phenoCol Y1"
        def ped = pheno_files.find{it.name.endsWith(".plink.ped")}
"""

mkdir -p TMP/OUT
set -x


regenie \\
  --step 1 \\
  --pgen \$(basename ${pgen} .pgen) \\
  --phenoFile ${ped} \\
  --phenoColList `head -n 1 ${ped} | cut -f4- |tr "\t" ","` \\
  --covarFile "${covariates}" \\
  --extract '${keep_rs}' \\
  ${args} \\
  --lowmem \\
  --lowmem-prefix TMP/regenie_tmp_preds \\
  --threads ${task.cpus} \\
  --out "${prefix}"

find ./
{meta1.id}step1

cat << EOF > versions.yml
${task.process}:
    regenie: TODO
EOF
"""
stub:
    def prefix = task.ext.prefix?:"{meta1.id}step1"
"""
touch versions.yml ${prefix}_pred.list ${prefix}.log
"""
}

