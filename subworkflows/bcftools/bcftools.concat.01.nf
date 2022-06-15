/**
*/

workflow BCFTOOLS_CONCAT {
	take:
		meta
		vcfs
		bed
	main:
		vers_ch = concat_version(meta, vcfs, bed)
		each_list_ch = sqrt_files(meta,vcfs).splitCsv()
		concat0_ch = concat0(meta,each_list_ch,bed)
		concat1_ch = concat1(meta,concat0_ch.vcf.collect())
	emit:
		vcf = concat1_ch.vcf
		version = vers_ch.version
	}


process concat0 {
tag "${file(vcfs).name}"
input:
	val(meta)
	val(vcfs)
	val(bed)
output:
	path("concat.0.bcf"),emit:vcf	
	path("concat.0.bcf.csi"),emit:csi
script:
	"""
	hostname 1>&2
	module load ${getModules("bcftools")}

	bcftools concat --threads ${task.cpus} \
		${bed.isEmpty()?"":"--regions-file \"${bed}\""} |\
		--no-version --allow-overlaps --remove-duplicates \
		-O b -o "concat.0.bcf" --file-list "${vcfs}"

	bcftools index --threads ${task.cpus}  "concat.0.bcf"
	"""
	}

process concat_version {
	executor "local"
	input:
		meta
		vcfs
		bed
	output:
		path("version.xml"),emit:version
	script:
	"""
	cat << EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">concat vcf(s) using bcftools</entry>
		<entry key="vcfs">${vcfs}</entry>
		<entry key="bed">${bed}</entry>
	</properties>
	EOF
	"""
	}

process concat1 {
tag "N=${L.size()}"
input:
	val(meta)
	val(L)
output:
	path("${prefix}concat${suffix}"),emit:vcf	
script:
	prefix = getKeyValue(meta,"prefix","")
	sufffix = getKeyValue(meta,"suffix",".bcf")

	if(L.size()==1)
	"""
	hostname 1>&2
	module load ${getModules("bcftools")}
	
	bcftools sort -T . \
		-O "${suffix.contains("b")?"b":"z"} \
		-o "${prefix}concat${suffix.contains("b")?".bcf":".vcf.gz"}" \
		"${L[0]}"

	bcftools index --threads ${task.cpus} \
		${prefix}concat${suffix.contains("b")?"":"--tbi"} \
		"${prefix}concat${suffix.contains("b")?".bcf":".vcf.gz"}"

	"""
	else
	"""
	hostname 1>&2
	module load ${getModules("bcftools")}
	set -o pipefail

cat << EOF > jeter.list
${L.join(" ")}
EOF

	bcftools concat --threads ${task.cpus} \
		--no-version --allow-overlaps --remove-duplicates \
		-O u --file-list jeter.list |\
		bcftools sort -T . \
			-O "${suffix.contains("b")?"b":"z"} \
			-o "${prefix}concat${suffix.contains("b")?".bcf":".vcf.gz"}" 

	bcftools index --threads ${task.cpus} \
		${prefix}concat${suffix.contains("b")?"":"--tbi"} \
		"${prefix}concat${suffix.contains("b")?".bcf":".vcf.gz"}"

	rm jeter.list
	"""
	}

