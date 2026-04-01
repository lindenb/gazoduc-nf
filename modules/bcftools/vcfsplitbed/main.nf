

process VCF_SPLIT_BED
tag "${meta.id}"
label "process_single"
input:
	tuple val(meta1),path(bed)
	tuple val(meta ),path(vcf),path(tbi)
output:
	tuple val(meta),path("*.vcf.gz",arity:"0..*"),emit:vcf
script:
	def prefix = task.ext.prefix?:"${meta.id}.${meta1.id}"
	def args1 = task.ext.args1?:""
"""
mkdir -p TMP

${bed.name.endsWith(".gz")?"gunzip -c":"cat"} "${bed}" |\\
	grep -v '^(#|browser|track)' |\\
	awk -F '\t' '{N=\$4;if(NF< 4 || N=="") N=sprintf("%s_%s_%s",\$1,\$2,\$3); printf("%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,N);  }' |\\
	sort -T TMP -t '\t' -k4,4 |\\
	awk -F '\t' 'BEGIN {PREV="";OFS="\t";} {if(PREV!=\$4) {if(NR>1)printf("\\\n"); PREV=\$4; printf("%s",PREV); } printf("|%s,%s,%s",\$1,\$2,\$3);} END {if(PREV!="") {printf("\\n");} } ' > TMP/by_gene.txt

cat TMP/by_gene.txt | while read L
do
	echo "${L}" > TMP/jeter.txt
	NAME=`cut -d '|' TMP/jeter.txt`
	test -n "\${NAME}"
	cut -d '|' -f 2- TMP/jeter.txt | tr "|" "\\n" | tr "," "\t" |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\\
		bedtools merge > TMP/jeter.bed

	test -s TMP/jeter.bed

	bcftools view ${args1}  --regions-file TMP/jeter.bed -O z -o TMP/jeter.vcf.gz

	if test `bcftools query -f '.' TMP/jeter.vcf.gz | wc -c` -gt 0 
	then
		mv TMP/jeter.vcf.gz TMP/${prefix}.\${NAME}.vcf.gz
	else
		rm TMP/jeter.vcf.gz
	fi
done

find TMP -type f  -name "*.vcf.gz"  -exec mv -t "\${PWD}/" '{}' ';'

touch versions.yml
"""

stub:
"""
touch versions.yml
"""
}
