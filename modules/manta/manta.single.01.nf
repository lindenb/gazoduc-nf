process MANTA_SINGLE_01 {
    tag "${name} ${file(bam).name}"
    afterScript "rm -rf TMP"
    cache 'lenient'
    errorStrategy 'finish'
    cpus 16
    input:
	val(meta)
	val(reference)
	tuple val(name),val(bam)
    output:
    	tuple val(name),val(bame),path("${name}.txt"),emit:output
	path("version.xml"),emit:version
    script:
	def prefix = meta.getOrDefault("prefix","")
	"""
	hostname 1>&2
	module load ${getModules("manta")}
	mkdir -p TMP

	configManta.py  --bam "${bam}" --referenceFasta "${reference}" \
		--runDir "TMP"

	
	./TMP/runWorkflow.py --quiet -m local -j ${task.cpus}
	
	rm -rf ./TMP/workspace

	# change name to sample
	find ./TMP -type f -name "*.vcf.gz" \
		-printf 'mv -v %p  ${prefix}${name}.%f\\n mv -v %p.tbi ${prefix}${name}.%f.tbi\\n' |bash 
	
	ls *.vcf.gz
	find \${PWD}  -maxdepth 1  -type f -name "*.vcf.gz" -o -name "*.vcf.gz.tbi" |\
		awk '{printf("${name}\t${bam}\t%s\\n");}' > ${name}.txt
	grep -m1 vcf ${name}.txt

#################################################################################################

cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">${description}</entry>
	<entry key="date">\$(date)</entry>
	<entry key="sample">${name}</entry>
	<entry key="bam">${bam}</entry>
	<entry key="manta.version">\$(configManta.py --version)</entry>
</properties>
EOF


	"""
	}
