/*

Copyright (c) 2022 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
The MIT License (MIT)
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
include {moduleLoad;getVersionCmd} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'

workflow MANTA_MULTI_SV01 {
	take:
		meta
		reference
		bams
		bed
	main:
		version_ch = Channel.empty()

		config_ch = MANTA_CONFIG(meta,reference,bams,bed)
		version_ch= version_ch.mix(config_ch.version)

		
		manta_ch = MANTA_MULTI(meta,reference,config_ch.output)
		version_ch= version_ch.mix(manta_ch.version)
	
		
		version_ch = MERGE_VERSION(meta,"manta","manta",version_ch.collect())
		html = VERSION_TO_HTML(params, version_ch.version)

		zip_ch = ZIPIT(meta, version_ch, html.html, manta_ch.output)	
	emit:
		version = version_ch
		zip = zip_ch.zip
		//vcf = manta_ch.vcf
		//index = manta_ch.index
	}

process MANTA_CONFIG {
    tag "${bams.name}"
    input:
	val(meta)
	val(reference)
	path(bams)
	path(bed)

    output:
    	path("path.txt"),emit:output
	path("version.xml"),emit:version

    script:
	def workingdir = "${workflow.workDir}/${meta.prefix?:""}MANTA"
	log.info("MANTA WORKDIR ${workingdir}")
	"""
	hostname 1>&2
	${moduleLoad("manta htslib")}
	

	rm -rvf "${workingdir}"

	if ${bed.name.equals("NO_FILE")} ; then
		awk -F '\t' '(\$1 ~ /^(chr)?[0-9XY]+\$/ ) {printf("%s\t0\t%s\\n",\$1,\$2)}' '${reference}.fai' |\
		sort -T . -t '\t' -k1,1 -k2,2n > manta.bed
	else
		${bed.name.endsWith(".gz")?"gunzip -c":"cat"} ${bed} |\
		sort -T . -t '\t' -k1,1 -k2,2n > manta.bed
	
	fi

	bgzip manta.bed
	tabix -p bed manta.bed.gz

	configManta.py  `awk '{printf("--bam %s ",\$0);}' "${bams}" ` \
		--referenceFasta "${reference}" \
		--callRegions manta.bed.gz \
		--runDir "${workingdir}"

	echo "${workingdir}" > path.txt

#########################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">Manta configuration</entry>
	<entry key="bams">${bams}</entry>
	<entry key="workdir">${workingdir}</entry>
	<entry key="manta.version">\$(configManta.py --version)</entry>
	<entry key="versions">${getVersionCmd("tabix")}</entry>
</properties>
EOF
	"""
	}

process MANTA_MULTI {
    cache 'lenient'
    errorStrategy 'retry'
    clusterOptions = "-S /bin/bash -q max-7d.q "
    maxRetries 20
    cpus 16
    memory "20g"
    input:
	val(meta)
	val(reference)
	path(workingdir)
    output:
	path("paths.txt"),emit:output
	path("version.xml"),emit:version
    script:
	def prefix = meta.prefix?:""
	"""
	hostname 1>&2
	${moduleLoad("manta bcftools")}
		
	WD=`cat "${workingdir}"`
	test -d "\${WD}"

	"\${WD}/runWorkflow.py" --quiet -m local -j ${task.cpus}
	

	for X in diploidSV candidateSmallIndels candidateSV
	do
		bcftools view -O b -o "${meta.prefix?:""}\${X}.bcf" "\${WD}/results/variants/\${X}.vcf.gz"
		bcftools index "${meta.prefix?:""}\${X}.bcf"
		echo "\${PWD}/${meta.prefix?:""}\${X}.bcf" >> paths.txt
		echo "\${PWD}/${meta.prefix?:""}\${X}.bcf.csi" >> paths.txt
	done


#################################################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="workdir">\${WD}</entry>
	<entry key="versions">${getVersionCmd("bcftools")}</entry>
	<entry key="manta.version">\$(configManta.py --version)</entry>
</properties>
EOF
	"""
	}


process ZIPIT {
input:
	val(meta)
	path(version)
	path(html)
	path(paths)
output:
	path("${meta.prefix?:""}manta.zip"),emit:zip
script:
"""
cat "${paths}" > x.list
echo "${version.toRealPath()}" >> x.list
echo "${html.toRealPath()}" >> x.list

zip -9 -@ -j "${meta.prefix?:""}manta.zip" < x.list

rm x.list
"""
}
