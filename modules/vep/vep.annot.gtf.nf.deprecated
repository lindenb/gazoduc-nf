
process VEP_ANNOT_GTF_01 {
tag "${file(vcf).name}"
afterScript "rm -f jeter.bcf"
memory "5g"
input:
	val(meta)
	val(vcf)
	val(bed)
	val(gff)
output:
	path("vep.bcf"),emit:vcf
	path("vep.bcf.csi"),emit:csi
	path("version.xml"),emit:version
script:
	def preBcftoolExtraArgs = getKeyValue(meta,"pre_bcftools_args","")
	def postBcftoolExtraArgs = getKeyValue(meta,"post_bcftools_args","")
"""
hostname 1>&2
set -o pipefail
module load ${getModules("vep")}

bcftools view ${preBcftoolExtraArgs} ${bed.isEmpty()?"":"--regions-file \""+bed+"\""} "${vcf}" |\
	vep --verbose --format vcf --force_overwrite --gtf "${gtf}"  --fasta "${reference}" --use_given_ref --vcf |\
	bcftools view ${postBcftoolExtraArgs} -o jeter.bcf -O b

bcftools index jeter.bcf

mv jeter.bcf "vep.bcf"
mv jeter.bcf.csi "vep.bcf.csi"


##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">annot with VEP</entry>
	<entry key="vcf">${vcf}</entry>
	<entry key="bed">${bed}</entry>
</properties>
EOF
"""
stub:
"""
touch vep.bcf vep.bcf.csi
echo "<properties/>" > version.xml
"""
}
