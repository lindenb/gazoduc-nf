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

process REGENIE_STEP2 {
label "process_single"
array 100
tag "${meta.id}"
conda "${moduleDir}/../../conda/regenie.yml"
afterScript "rm -rf TMP"
input:
        tuple val(meta1 ),path(pgen),path(psam),path(pvar)
        tuple val(meta2),path(covariates)
        tuple val(meta3),path(plink_ped)
	tuple val(meta4),path(loco)
	tuple val(meta ),path(annot),path(setfile),path(mask),path(aff)
output:
        tuple val(meta),path("*.regenie.gz"),emit:output
	tuple val(meta),path("*masks.snplist.gz"),emit:masks_snplist
script:
	 def prefix = task.ext.prefix?:"${meta.id}step2"
"""

mkdir -p TMP/OUT
set -x

gunzip -c "${annot}" > TMP/annot.txt
gunzip -c "${setfile}" > TMP/setfile.txt
gunzip -c "${mask}" > TMP/mask.txt
gunzip -c "${aff}" > TMP/aaf.txt

regenie \\
  --step 2 \\
  --pgen ${pgen.baseName} \\
  --phenoFile ${plink_ped} \\
  --covarFile "${covariates}" \\
  --pred ${loco} \\
  --mask-def TMP/mask.txt \\
  --set-list TMP/setfile.txt \\
  --anno-file TMP/annot.txt \\
  ${params.use_aaf_file?"--aaf-file TMP/aaf.txt":""} \\
  --phenoCol ${params.status} \\
  --bt \\
  --bsize 1000 \\
  --lowmem \\
  --lowmem-prefix TMP/regenie_tmp_preds \\
  --threads ${task.cpus} \\
  --out "step2.${contig}" \\
  --aaf-bins ${params.freq} \\
  --vc-maxAAF ${params.vc_maxAAF} \\
  --bsize 200 \\
  --vc-tests "${params.vc_tests}" \\
  --check-burden-files \\
  --weights-col ${params.weight_column?:4} \\
  --write-mask-snplist \\
  --firth --approx \\
  --pThresh 0.01

gzip --best *masks.snplist
gzip --best *.regenie

cat << EOF > versions.yml
${task.process}:
    regenie: TODO
EOF
"""
stub:
    def prefix = task.ext.prefix?:"${meta.id}step2"
"""
touch versions.yml 
"""
}

