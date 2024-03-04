include {slurpJsonFile;moduleLoad} from '../../modules/utils/functions.nf'


def TAG="FIRST"
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
	path("OUTPUT/${TAG}.json"),emit:output
	path("OUTPUT/${TAG}.count"),emit:count
script:
	if(!interval.isEmpty() && !bed.name.equals("NO_FILE")) throw new IllegalArgumentException("both bed and interval are defined.");
	def args = " --regions-file \"TMP/FIRST.bed\""
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP OUTPUT

if ${interval.isEmpty()} ; then
	ln -s "${bed}" TMP/FIRST.bed
else
	echo "${interval}" | awk -F '[:-]' '{printf("%s\t%d\t%s\\n",\$1,int(\$2)-1,\$3);}' > TMP/${TAG}.bed
fi

bcftools view ${args} -O b -o TMP/${TAG}.bcf "${vcf}"
bcftools index --force TMP/${TAG}.bcf

cat << EOF > TMP/${TAG}.json
{
"vcf"   : "\${PWD}/OUTPUT/${TAG}.bcf",
"index" : "\${PWD}/OUTPUT/${TAG}.bcf.csi",
"bed"   : "\${PWD}/OUTPUT/${TAG}.bed"
}
EOF


### count variants
bcftools query -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/${TAG}.count
mv -v TMP/${TAG}.* OUTPUT
"""
}
