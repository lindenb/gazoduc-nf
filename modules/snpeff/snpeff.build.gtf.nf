include {getKeyValue;getModules} from '../utils/functions.nf'

/**
 * build a SNPEFF database from a GTF file
 */
process SNPEFF_BUILD_GTF {
tag "${file(gtf).name} ${file(fasta).name}"
afterScript "rm -f snpeff.errors"
memory "10g"
input:
    val(meta)
    val(fasta)
    val(gtf)
output:
	tuple val("${meta.dbName?:file(gtf).getSimpleName()}"),path("snpEff.config"),emit:snpeffdb
	path("version.xml"),emit:version
script:
	def dbName = meta.dbName?:file(gtf).getSimpleName()
"""
	hostname 1>&2
	module load ${getModules("snpEff jvarkit")}
	set -o pipefail
	
	mkdir -p "data/${dbName}"
	ln -s "${fasta}" "data/${dbName}/sequences.fa"

	# convert gtf chromosomes
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar \${JVARKIT_DIST}/bedrenamechr.jar \
		-R "${fasta}" --column 1 --convert SKIP  "${gtf}" > "data/${dbName}/genes.gtf"
	test -s "data/${dbName}/genes.gtf"

	# write snpEff contig
	cat <<- EOF > snpEff.config
	data.dir = \${PWD}/data/
	${dbName}.genome = Human
	${dbName}.reference = ${fasta}
	EOF

	# build database
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar \${SNPEFF_JAR}  build -gtf22 -v "${dbName}"

	test -s "data/${dbName}/snpEffectPredictor.bin"

	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">build custom snpeff from gtf</entry>
		<entry key="gtf">${gtf}</entry>
		<entry key="fasta">${fasta}</entry>
		<entry key="output">\${PWD}/snpEff.config</entry>
		<entry key="bedrenamechr">\$(java  -jar \${JVARKIT_DIST}/bedrenamechr.jar --version)</entry>
		<entry key="snpeff">\$(java -jar \${SNPEFF_JAR} -version)</entry>
	</properties>
	EOF
"""

stub:
"""
touch "${meta.dbName?:file(gtf).getSimpleName()}" "snpEff.config"
echo "<properties/>" > version.xml
"""
}
