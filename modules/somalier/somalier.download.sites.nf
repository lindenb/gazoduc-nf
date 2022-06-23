
process SOMALIER_DOWNLOAD_SITES {
tag "${file(reference).name}"
input:
	val(meta)
	val(reference)
output:
	path("sites.vcf.gz"), emit: vcf
	path("sites.vcf.gz.tbi")
	path("version.xml"), emit:version
script:
	def url=(isHg19(reference)?"https://github.com/brentp/somalier/files/3412453/sites.hg19.vcf.gz":
		(isHg38(reference)?"https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz":
		""))
"""
hostname 1>&2
${moduleLoad("jvarkit bcftools")}
set -o pipefail
set -x

wget -O - "${url}" |\
	gunzip -c |\
	java -jar \${JVARKIT_DIST}/vcfsetdict.jar -R "${reference}"  --onNotFound SKIP |\
	bcftools sort -T . -o sites.vcf.gz -O z

bcftools index -t sites.vcf.gz

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">download VCF sites for somalier</entry>
        <entry key="url">${url}</entry>
	<entry key="bcftools.version">\$( bcftools --version-only)</entry>
</properties>
EOF
"""
}
