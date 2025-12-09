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
process MERGE_BEDS {
label "process_single"
tag "${meta.id}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/goleft.yml"

input:
    tuple val(meta),path("BED/*")
output:
    tuple val(meta),path("*.bed.gz"),path("*.bed.gz.tbi"),emit:bed
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:meta.id+".merge"
"""
hostname 1>&2
mkdir -p TMP
set +o pipefail

find BED/ -name "*.bed.gz" | while read F
do
	if ! test -f TMP/jeter.bed
	then

		gunzip -c "\${F}" | head -n 1 > TMP/header.txt
		gunzip -c "\${F}" | tail -n +2 | cut -f 1,2,3 > TMP/signature1.txt
		gunzip -c "\${F}" | tail -n +2 > TMP/jeter.bed

	else
		gunzip -c "\${F}" | tail -n +2 | cut -f 1,2,3 > TMP/signature2.txt
		# check all files have the same intervals in the same order
		cmp TMP/signature1.txt TMP/signature2.txt

		paste TMP/header.txt <(gunzip -c "\${F}" | head -n 1 | cut -f4-)  > TMP/jeter2.txt
		mv TMP/jeter2.txt TMP/header.txt

		paste TMP/jeter.bed <(gunzip -c "\${F}" | tail -n +2 | cut -f4-) > TMP/jeter2.txt
		mv TMP/jeter2.txt TMP/jeter.bed
	fi
done

cat TMP/header.txt TMP/jeter.bed > TMP/jeter2.txt
mv TMP/jeter2.txt TMP/jeter.bed

# check there is only one number of cols
test \$(awk '{print NF}' TMP/jeter.bed | uniq | sort | uniq | wc -l) -eq 1

mv TMP/jeter.bed "TMP/${prefix}.bed"


bgzip "TMP/${prefix}.bed"

tabix -f -p bed "TMP/${prefix}.bed.gz"

mv TMP/*.gz ./
mv TMP/*.gz.tbi ./

touch versions.yml
"""

stub:
 def prefix = task.ext.prefix?:meta.id+".merge"
"""
touch versions.yml  ${prefix}.bed.gz.tbi
echo -e "chrom\tstart\tend\tS1\tS2\tS3\tS4" > ${prefix}.bed
gzip ${prefix}.bed
"""
}
