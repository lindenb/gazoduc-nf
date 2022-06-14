
process vcf2bed {
tag "${file(vcfs).name}"
input:
	val(meta)
	val(vcf)
output:
	path("vcf2bed.bed"),emit:bed
	path("version.xml"),emit:version
script:

	if(vcf.endsWith(".list"))
	"""
	hostname 1>&2
	set -o pipefail
	module load ${getModules("bcftools")}

	cat "${vcf}" | while read V
	do
		bcftools index -s "\${V}" | awk -F '\t' -vV=\${V} '{printf("%s\t0\t%s\t%s\\n",\$1,\$2,V);}'
	done | sort -T . -t '\t' -k1,1 -k2,2n | uniq > vcf2bed.bed


	
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">extract contigs in multiple VCFs using bcftools</entry>
		<entry key="input">${vcf}</entry>
	</properties>
	EOF
	"""
	else
	"""
	hostname 1>&2
	set -o pipefail
	module load ${getModules("bcftools")}

	bcftools index -s "${vcf}" |\
		awk -F '\t' '{printf("%s\t0\t%s\t${vcf}\\n",\$1,\$2);}' |\
		sort -T . -t '\t' -k1,1 -k2,2n |\
		uniq > vcf2bed.bed


	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">extract contigs in one VCFs using bcftools</entry>
		<entry key="input">${vcf}</entry>
	</properties>
	EOF
	"""
	}
