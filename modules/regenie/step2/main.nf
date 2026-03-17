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
include  { isBlank } from '../../../modules/utils/functions.nf'
include  { verify  } from '../../../modules/utils/functions.nf'


process REGENIE_STEP2 {
label "process_single"
label "array100"
tag "${meta.id}"
conda "${moduleDir}/../../../conda/regenie.yml"
afterScript "rm -rf TMP"
input:
    tuple val(meta1 ),path(pgen),path(pvar),path(psam)
    tuple val(meta2),path(covariates)
    tuple val(meta3),path(plink_ped)
    tuple val(meta4),path(loco)
    tuple val(meta ),path(annot),path(setfile),path(aff),path(mask)
output:
    tuple val(meta),path("*.regenie.gz"),emit:regenie
    tuple val(meta),path("*masks.snplist.gz"),emit:masks_snplist
    path("versions.yml"),emit:versions
script:     
   def prefix = task.ext.prefix?:"${meta.id}.step1"
   def input_arg = (pgen.name.endsWith(".pgen")?"--pgen":(pgen.name.endsWith(".bgen")?"--bgen":"--bed"))
   def use_aaf_file = task.ext.use_aaf_file?:true
   def phenoCol = task.ext.phenoCol?:"status"
   verify(!isBlank(task.ext.aaf_bins),"${task.process} needs task.ext.aaf_bins")
   def aaf_bins = task.ext.aaf_bins?:"-1"
   verify(!isBlank(task.ext.vc_maxAAF),"${task.vc_maxAAF} needs task.ext.vc_maxAAF")
   def vc_maxAAF = task.ext.vc_maxAAF?:"-1"
   verify(!isBlank(task.ext.vc_tests),"${task.vc_tests} needs task.ext.vc_tests")
   def vc_tests = task.ext.vc_tests?:"-1"
   def weight_column = task.ext.weight_column?:"4"
"""

mkdir -p TMP/OUT
set -x

echo '${phenoCol} ${loco}' > TMP/regenie_step1_pred.list
gunzip -c "${annot}" > TMP/annot.txt
gunzip -c "${setfile}" > TMP/setfile.txt
gunzip -c "${mask}" > TMP/mask.txt
gunzip -c "${aff}" > TMP/aaf.txt

regenie \\
  --step 2 \\
  ${input_arg} ${pgen.baseName} \\
  --phenoFile ${plink_ped} \\
  --covarFile "${covariates}" \\
  --pred TMP/regenie_step1_pred.list \\
  --mask-def TMP/mask.txt \\
  --set-list TMP/setfile.txt \\
  --anno-file TMP/annot.txt \\
  ${use_aaf_file?"--aaf-file TMP/aaf.txt":""} \\
  --phenoCol ${phenoCol} \\
  --bt \\
  --bsize 1000 \\
  --lowmem \\
  --lowmem-prefix TMP/regenie_tmp_preds \\
  --threads ${task.cpus} \\
  --out "${prefix}" \\
  --aaf-bins ${aaf_bins} \\
  --vc-maxAAF ${vc_maxAAF} \\
  --bsize 200 \\
  --vc-tests "${vc_tests}" \\
  --check-burden-files \\
  --weights-col ${weight_column} \\
  --write-mask-snplist \\
  --firth --approx \\
  --pThresh 0.01

gzip --best *masks.snplist
gzip --best *.regenie

cat << EOF > versions.yml
${task.process}:
    regenie: \$(regenie --help |grep REGENIE -m1 | awk '{print \$3}')
EOF
"""
stub:
    def prefix = task.ext.prefix?:"${meta.id}step2"
"""
touch versions.yml 
"""
}

