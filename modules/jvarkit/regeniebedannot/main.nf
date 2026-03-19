
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
include { verify   } from '../../../modules/utils/functions.nf'
include { isBlank  } from '../../../modules/utils/functions.nf'


process REGENIE_BED_ANNOT {
tag "${meta.id}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"   
afterScript "rm -rf TMP"
input:
    tuple val(meta1),path(select_bed)
    tuple val(meta ),path(vcf),path(vcf_tbi)
output:
    tuple val(meta),path("*.tsv.gz"),emit:tsv
	path("versions.yml"),emit:versions
script:
    def args1 = task.ext.args1?:""
    def args3 = task.ext.args3?:""
    def awk_expr = task.ext.awk_expr?:"(1==1)"
	def jvm = task.ext.jvm?:"-Djava.io.tmpdir=TMP "
	def jvarkit = task.ext.jvarkit?:"java -jar  ${jvm} \${HOME}/jvarkit.jar"
	def freq = task.ext.freq?:""
    verify(!isBlank("${freq}"),"${task.process} freq is blank ?")
    def prefix = task.ext.prefix?:"${meta.id}.${meta1.id}"
	def min_length= task.ext.min_bed_length?:"0" //(params.min_bed_length?:0)
"""
mkdir -p TMP
set -x

# extract DICT only
bcftools view --header-only  -O z -o TMP/dict.vcf.gz '${vcf}'

# rename contig in user bed
${select_bed.name.endsWith(".gz")?"gunzip -c ":"cat"} "${select_bed}" |\\
	${jvarkit} bedrenamechr -R TMP/dict.vcf.gz --column 1 |\\
    awk -F '\t' '(${awk_expr})' > TMP/jeter.bed

# prevent empty file
if test ! -s TMP/jeter.bed
then
	echo -e "XXXXX\t0\t1" > TMP/jeter.bed
fi

# give a chance to select a whole chromosome with args1
bcftools view ${args1} -O u '${vcf}' |\\
   bcftools view -O v --targets-file TMP/jeter.bed --targets-overlap 2 |\\
	${jvarkit} regeniebedannot \\
        --bed "${select_bed}" \\
		--min-length '${min_length}' \\
        ${args3} \\
		-f "${freq}" |\\
        gzip --best > TMP/jeter.tsv.gz

mv -v TMP/jeter.tsv.gz ${prefix}.tsv.gz
    
cat << EOF > versions.yml
${task.process}:
	jvarkit: "\$(${jvarkit} --version)"
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
EOF
"""

stub:
	def freq = task.ext.freq?:""
    verify(!isBlank("${freq}"),"${task.process} freq is blank ?")
    def prefix = task.ext.prefix?:"${meta.id}.${meta1.id}"
"""

touch versions.yml ${prefix}.tsv.gz
"""
}
