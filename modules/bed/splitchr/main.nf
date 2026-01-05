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

process SPLIT_BED_PER_CHROMOSOME {
tag "${meta.id} ${bed.name}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta),path(bed)
output:
    tuple val(meta),path("*.{bed,bed.gz}",arity:"0..*"),optional:true,emit:beds
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"${meta.id}.percontig"
    def with_bgzip = task.ext.with_bgzip?:false
"""
mkdir -p TMP

${bed.endsWith(".gz")?"gunzip -c":"cat"} "${bed}" |\\
	grep -vE '^(browser|track|#)' |\\
	sort  -S ${task.memory.kilo} -t '\t' -T TMP -k1,1 -k2,2n |\\
	awk -F '\t' '
	BEGIN {
		PREV = "";
		OFS = "\t";
		}
	{
	C=\$1;
	if(PREV != C) {
		if(PREV!="") close(f);
		PREV=C;
		f = sprintf("TMP/${prefix}.%s.bed",C);
		}
	print >> f;
	}'
	

if ${with_bgzip}
then
	find ./TMP -name "*.bed" -exec bgzip '{}' ';'
	mv TMP/*.bed.gz ./ || true
else
	mv TMP/*.bed ./ || true
fi


touch versions.yml
"""
stub:
	def prefix = task.ext.prefix?:"${meta.id}"
"""
for F in 1 2 3 4 5 6
do
	echo "chr\t\${F}\t0\t100000" > "${prefix}.\${F}.bed"
done
touch versions.yml
"""
}
