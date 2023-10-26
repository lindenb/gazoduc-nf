/*

Copyright (c) 2023 Pierre Lindenbaum

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
nextflow.enable.dsl=2

def gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults().putReference();

gazoduc.make("cases","NO_FILE").
        description("file containing the path to multiple bam files for the cases").
        required().
        existingFile().
        put()

gazoduc.make("controls","NO_FILE").
        description("file containing the path to multiple bam files for the controls").
        required().
        existingFile().
        put()

gazoduc.putCondaEnv()


include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {parseBoolean;isHg19;getVersionCmd;moduleLoad;runOnComplete} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {SAMTOOLS_CASES_CONTROLS_01} from '../../subworkflows/samtools/samtools.cases.controls.01.nf'
include {SIMPLE_ZIP_01} from '../../modules/utils/zip.simple.01.nf'

if( params.help ) {
    gazoduc.usage().
		name("expansion hunter de novo").
		description("expansion hunter de novo").
		print();
    exit 0
    }
else
	{
	gazoduc.validate()
	}



workflow {
	ch1 = EXPANSION_HUNTER_DE_NOVO_01(params,params.reference,file(params.cases),file(params.controls))
	html = VERSION_TO_HTML(params,ch1.version)
	to_zip = Channel.empty().
		mix(ch1.version).
		mix(ch1.output).
		mix(html.html)
	SIMPLE_ZIP_01([:],to_zip.collect())
	}

runOnComplete(workflow);

workflow EXPANSION_HUNTER_DE_NOVO_01 {
	take:
		meta
		reference
		cases
		controls
	main:

		version_ch = Channel.empty()

	
		cases_controls_ch = SAMTOOLS_CASES_CONTROLS_01([:],reference,cases,controls)
		version_ch= version_ch.mix(cases_controls_ch.version)


		xhunter = DOWNLOAD_XHUNTER(meta)
		version_ch = version_ch.mix(xhunter.version)

		each_case_control_ch = cases_controls_ch.output.splitCsv(header:false,sep:'\t')

		runxhunter_ch = APPLY_XHUNTER(meta, reference, xhunter.xhunter, each_case_control_ch)
		version_ch = version_ch.mix(runxhunter_ch.version)

		merge_ch= XHUNTER_MERGE(meta, reference, xhunter.xhunter,runxhunter_ch.output.collect())
		version_ch = version_ch.mix(merge_ch.version)

		methods = Channel.of("motif","locus")

		xcc_ch = XHUNTER_CASE_CONTROL(meta,reference, xhunter.xhunter, merge_ch.json, merge_ch.manifest, methods)
		version_ch = version_ch.mix(xcc_ch.version)

		version_ch = MERGE_VERSION(meta, "ExpansionHunterDeNovo", "Expansion hunter de-novo", version_ch.collect())
	emit:
		version = version_ch
		output = xcc_ch.output
	}

process DOWNLOAD_XHUNTER {
	executor "local"
	afterScript "rm -f jeter.tar.gz"
	input:
		val(meta)
	output:
		path("ExpansionHunterDenovo-linux_x86_64/bin/ExpansionHunterDenovo"),emit:xhunter
		path("version.xml"),emit:version
	script:
		def xversion = meta.expansion_hunter_version?:"0.9.0"
	"""
	hostname 1>&2
	set -x
	wget -O jeter.tar.gz "https://github.com/Illumina/ExpansionHunterDenovo/releases/download/v${xversion}/ExpansionHunterDenovo-v${xversion}-linux_x86_64.tar.gz"
	tar xvfz jeter.tar.gz
	chmod +x "ExpansionHunterDenovo-v${xversion}-linux_x86_64/bin/ExpansionHunterDenovo"
	# rename directory
	mv "ExpansionHunterDenovo-v${xversion}-linux_x86_64"  "ExpansionHunterDenovo-linux_x86_64"

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">Download expansion hunter</entry>
	<entry key="expansion.hunter.version">${xversion}</entry>
	<entry key="version">${getVersionCmd("wget")}</entry>
</properties>
EOF
	"""
	}

process APPLY_XHUNTER {
	tag "${sample}=${status} ${file(bam).name}"
	afterScript "rm -rf TMP"
	memory "3g"
	cpus 1
	input:
		val(meta)
		val(reference)
		path(executable)
		tuple val(sample),val(bam),val(status)
	output:
		path("manifest.txt"),emit:output
		path("version.xml"),emit:version
	script:

	"""
	hostname 1>&2
        mkdir -p TMP

	${executable.toRealPath()} profile \
		--reads "${bam}" \
		--reference "${reference}" \
		--output-prefix "TMP/jeter" 1>&2
	

	mv -v TMP/jeter.str_profile.json  "${sample}.str_profile_json"
	echo "${sample}\t${status}\t\${PWD}/${sample}.str_profile_json" > manifest.txt

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">run expansion hunter denovo</entry>
	<entry key="sample">${sample}</entry>
	<entry key="bam">${bam}</entry>
</properties>
EOF
	"""
	}


process XHUNTER_MERGE {
	tag "N=${L.size()}"
	afterScript "rm -rf TMP"
	cpus 10
	input:
		val(meta)
		val(reference)
		path(executable)
		val(L)
	output:
		path("${meta.prefix?:""}multisample_profile.json"),emit:json
		path("manifest.txt"),emit:manifest
		path("version.xml"),emit:version
	script:
	"""
	hostname 1>&2
	${moduleLoad("bcftools")}
	set -o pipefail

mkdir -p TMP
cat << EOF > TMP/jeter.list
${L.join("\n")}
EOF

xargs -a TMP/jeter.list -L 50 cat > TMP/manifest.txt



	${executable.toRealPath()} merge \
		--manifest TMP/manifest.txt \
		--reference "${reference}" \
		--output-prefix "TMP/jeter" 1>&2

	mv TMP/jeter.multisample_profile.json "${meta.prefix?:""}multisample_profile.json"
	mv TMP/manifest.txt ./

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">merge json</entry>
	<entry key="json.count">${L.size()}</entry>
</properties>
EOF
	"""
	}


process XHUNTER_CASE_CONTROL {
	tag "${method}"
	afterScript "rm -rf TMP"
	conda "${params.conda}/SPLICEAI" //TODO fix this
	input:
		val(meta)
		val(reference)
		path(executable)
		path(json)
		path(manifest)
		val(method)
	output:
		path("${params.prefix?:""}${method}.tsv"),emit:output
		path("version.xml"),emit:version
	script:
	"""
	hostname 1>&2
	set -o pipefail


	\$(dirname "${executable.toRealPath()}" )/../scripts/casecontrol.py  ${method} \
		--manifest '${manifest}' \
		--multisample-profile '${json}' \
		--output jeter.tsv 


	mv jeter.tsv "${params.prefix?:""}${method}.tsv"

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">case control</entry>
	<entry key="method">${method}</entry>
</properties>
EOF
	"""
	}




