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
label "memory_50G"
tag "${meta.id}"
conda "${moduleDir}/../../../conda/regenie.yml"
afterScript "rm -rf TMP"
input:
        tuple val(meta2),path(covariates)
        tuple val(meta3),path(plink_ped)
        tuple val(meta4),path(keep_rs) // txt file or .bim file
        tuple val(meta ),path(pgen),path(psam),path(pvar)
output:
        tuple val(meta ),path("*_pred.list"),emit:pred_list
        tuple val(meta ),path("*.loco", arity: '1..*'),emit:loco
        tuple val(meta ),path("*.log"),emit:log
        path("versions.yml"),emit:versions
script:
        def prefix = task.ext.prefix?:"${meta.id}.step1"
        def args = task.ext.args?:"--bsize 1000 --bt --phenoCol Y1"
        def phenoColList = task.ext.phenoColList?:"status"
        def n_markers = 1000000
"""
mkdir -p TMP
set -x

if ${keep_rs?true:false}
then
        cp ${keep_rs} TMP/keep.markers.txt

        # from plink Bim file ?
        if ${(keep_rs?true:false) && keep_rs.name.endsWith(".bim")}
        then
                cut -f 2 TMP/keep.markers.txt  > TMP/jeter.txt
                mv TMP/jeter.txt TMP/keep.markers.txt

                head TMP/keep.markers.txt 1>&2
        fi



        #
        # it is not recommened to use more than ${n_markers} variants in step 1 of regenie
        #
      
        # shuffle and extract 
        awk '{printf("%d,%s\\n",int(rand()*1000000),\$0);}' TMP/keep.markers.txt |\\
                sort  -S ${task.memory.kilo} -t, -T TMP -k1,1n   |\\
                cut -d, -f2- |\\
                head -n ${n_markers} > TMP/jeter.txt
        mv TMP/jeter.txt TMP/keep.markers.txt

else

        awk -F '\t' '/^#/ {next;} {printf("%d,%s\\n",int(rand()*1000000),\$3);}' '${pvar}' |\\
                sort  -S ${task.memory.kilo} -t, -T TMP -k1,1n  |\\
                cut -d, -f2- > TMP/jeter.txt
        
        head -n ${n_markers} TMP/jeter.txt >  TMP/keep.markers.txt
       
fi

head TMP/keep.markers.txt 1>&2
wc -l TMP/keep.markers.txt 1>&2



regenie \\
  --step 1 \\
  --pgen ${pgen.baseName} \\
  --phenoFile ${plink_ped} \\
  --phenoColList ${phenoColList} \\
  --covarFile "${covariates}" \\
  --extract TMP/keep.markers.txt \\
  ${args} \\
  --lowmem \\
  --lowmem-prefix TMP/regenie_tmp_preds \\
  --threads ${task.cpus} \\
  --out "${prefix}"

find ./ 1>&2

cat << EOF > versions.yml
${task.process}:
    regenie: TODO
EOF
"""
stub:
    def prefix = task.ext.prefix?:"${meta.id}step1"
"""
touch versions.yml ${prefix}_pred.list ${prefix}.log
"""
}