
process DOWNLOAD_SOMALIER {
input:
	val(meta)
output:
	path("somalier"),emit:executable
	path("version.xml"),emit:version
script:
	def sommalier_version = getKeyValue(meta,"somalier_version","0.2.15")
	def url = "https://github.com/brentp/somalier/releases/download/v${somalier_version}/somalier"
"""
hostname 1>&2
wget -O somalier "${url}"
chmod +x somalier

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">download and somalier</entry>
        <entry key="somalier;version">${somalier_version}</entry>
        <entry key="url">${url}</entry>
</properties>
EOF
"""
}
