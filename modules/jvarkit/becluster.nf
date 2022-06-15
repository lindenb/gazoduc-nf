process bedcluster {
	tag "${file(bed).name}"
	input:
		val(meta)
		val(bed)
		val(reference)
	output:
		path("clusters.list"),emit:clusters
		path("version.xml"),emit:version
	script:
		def method = meta.method?:""
		def by_chrom = getBoolean(meta,"by_chromosome")
	"""
	hostname 1>&2
	set -o pipefail
	module load ${getModules("jvarkit")}

	test -z "${method}"

	mkdir BEDS
	java  -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar \${JVARKIT_DIST}/bedcluster.jar \
		-R "${reference}" \
		${by_chrom?"--chromosome":""} \
		${method} \
		-o BEDS "${bed}"

	find \${PWD}/BEDS -type f -name "*.bed" > clusters.list
	test -s clusters.list

	#####################################################################################

	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Cluster BED file</entry>
		<entry key="reference">${reference}</entry>
		<entry key="method">${method}</entry>
	</properties>
	EOF
	"""
	}
