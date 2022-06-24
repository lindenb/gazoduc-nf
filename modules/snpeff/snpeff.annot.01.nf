process SNPEFF_ANNOT_01 {
tag "${file(vcf).name} ${dbName}"
afterScript "rm -f jeter.bcf"
memory "5g"
input:
	val(meta)
	val(vcf)
	val(bed)
	val(snpEffCfg)
	val(dbName)
output:
	path("snpeff.bcf"),emit:vcf
	path("snpeff.bcf.csi"),emit:csi
	path("version.xml"),emit:version
script:
	def preBcftoolExtraArgs = getKeyValue(meta,"pre_bcftools_args","")
	def postBcftoolExtraArgs = getKeyValue(meta,"post_bcftools_args","")
"""
hostname 1>&2
set -o pipefail
module load ${getModules("snpEff bcftools")}

bcftools view ${preBcftoolExtraArgs} ${bed.isEmpty()?"":"--regions-file \""+bed+"\""} "${vcf}" |\
	java -Djava.io.tmpdir=. -jar "$${SNPEFF_JAR}" eff \
		-config "${snpEffCfg}" -interval "jeter.bed" \
		-nodownload -noNextProt -noMotif -noInteraction -noLog -noStats -chr chr -i vcf -o vcf ${dbName} |\
	bcftools view ${postBcftoolExtraArgs} -o jeter.bcf -O b

bcftools index jeter.bcf

mv jeter.bcf "snpeff.bcf"
mv jeter.bcf.csi "snpeff.bcf.csi"


##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">annot with snpEff</entry>
	<entry key="vcf">${vcf}</entry>
	<entry key="bed">${bed}</entry>
</properties>
EOF
"""

stub:
"""
touch snpeff.bcf snpeff.bcf.csi
echo "<properties/>" > version.xml
"""
}
