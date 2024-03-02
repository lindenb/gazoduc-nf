include {moduleLoad} from '../../modules/utils/functions.nf'

/**

We just copy the VCF, all the remaining steps will DELETE the previous VCF
to avoid having a GROWING storage.

*/
workflow STEP_FIRST {
	take:
		rows /** tuple vcf,vcf_index,bed,interval */
	main:
		step_ch = COPY(rows)
	emit:
		output = step_ch.output
		count = step_ch.count
}


process COPY {
tag "${vcf.name}"
afterScript "rm -rf TMP"
input:
	tuple path(vcf),path(index),path(bed),val(interval)
output:
	tuple path("OUTPUT/FIRST.bcf"),path("OUTPUT/FIRST.bcf.csi"),path("OUTPUT/FIRST.bed"),emit:output
	path("OUTPUT/count.tsv"),emit:count
script:
	if(!interval.isEmpty() && !bed.name.equals("NO_FILE")) throw new IllegalArgumentException("both bed and interval are defined.");
	def args = " --regions-file \"TMP/FIRST.bed\""
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP

if ${interval.isEmpty()} ; then
	ln -s "${bed}" TMP/FIRST.bed
else
	echo "${interval}" | awk -F '[:-]' '{printf("%s\t%d\t%s\\n",\$1,int(\$2)-1,\$3);}' > TMP/FIRST.bed
fi

bcftools view ${args} -O b -o TMP/FIRST.bcf "${vcf}"
bcftools index --force TMP/FIRST.bcf

### count variants
bcftools query -f '.'  TMP/FIRST.bcf | wc -c | awk '{printf("FIRST\t%s\\n",\$1);}' > TMP/count.tsv
mv -v TMP OUTPUT
"""
}
