/**
 * build a SNPEFF database from a GTF file
 */
process SNPEFF_BUILD_GTF {
tag "${file(gtf).name} ${file(fasta).name}"
memory "5g"
input:
    val(meta)
    val(fasta)
    val(gtf)
output:
	tuple val(meta.dbName?:file(gtf).getSimpleName()),path("snpEff.config"),emit:snpeffdb
	path("version.xml"),emit:version
script:
	def dbName = meta.dbName?:file(gtf).getSimpleName()
"""
	hostname 1>&2
	module load ${getModules("snpeff")}
	set -o pipefail
	
	mkdir -p "data/${dbName}"
	ln -s "${fasta}" "data/${dbName}/sequences.fa"

	# convert gtf chromosomes
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar \${JVARKIT_DIST}/bedrenamechr.jar \
		-f "${fasta}" --column 1 --convert SKIP  "${gtf}" > "data/${dbName}/genes.gtf"
	test -s "data/${dbName}/genes.gtf"

	# write snpEff contig
	cat <<- EOF > snpEff.config
	data.dir = \${PWD}/data/
	EOF

	# build database
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar \${SNPEFF_JAR}  build -gtf22 -v "${dbName}"

	test -s data/${dbName}/snpEffectPredictor.bin

	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">build custom snpeff from gtf</entry>
		<entry key="gtf">${gtf}</entry>
		<entry key="fasta">${fasta}</entry>
		<entry key="output">\${PWD}/snpEff.config</entry>
	</properties>
	EOF
"""
}
