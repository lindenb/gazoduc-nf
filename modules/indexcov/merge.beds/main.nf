
process MERGE_BEDS {
label "process_single"
tag "${meta.id}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/goleft.yml"
when:
    task.ext.when == null || task.ext.when
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
}
