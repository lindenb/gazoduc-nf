include {moduleLoad} from '../../modules/utils/functions.nf'
include {MULTIQC_01} from '../../modules/multiqc/multiqc.01.nf'
include {PARAMS_MULTIQC} from '../../modules/utils/params.multiqc.nf'

boolean isConfigFile(f) {
	def s =  f.name.toLowerCase();
	if(!(s.endsWith(".yaml") || s.endsWith(".yml"))) return false;
	return s.contains("config");
	}

workflow MULTIQC {
	take:
		files
	main:
                params4multiqc_ch = PARAMS_MULTIQC([:])

                mqc_ch = APPLY_MULTIQC(files.flatten().mix(params4multiqc_ch.output).collect())
	emit:
		zip = mqc_ch.zip
	}


process APPLY_MULTIQC {
	tag "N=${files.size()}"
	afterScript "rm -rf TMP"
	input:
		val(files)
	output:
		path("${params.prefix?:""}multiqc.zip"),emit:zip
		path("${params.prefix?:""}multiqc/${params.prefix?:""}multiqc_report_data"),emit:datadir
		path("version.xml"),emit:version
	script:
		def prefix = params.prefix?:""
		def extra = ""//meta.extra?:""
		def title = ""//meta.title?"--title \"${meta.title}\"":""
		def comment = ""//meta.comment?"--comment \"${meta.comment}\"":""
		def configs = files.findAll{isConfigFile(it)}.collect{"--config ${it}"}.join(" ")
	"""
		hostname 1>&2
		${moduleLoad("multiqc")}
		mkdir -p TMP

		export TMPDIR=\${PWD}/TMP

cat << EOF > TMP/jeter.list
${files.findAll{!isConfigFile(it)}.join("\n")}
EOF

		mkdir -p "${prefix}multiqc"

		export LC_ALL=en_US.utf8
		multiqc  --filename  "${prefix}multiqc_report.html" --no-ansi \
			${title}  \
			${comment}  \
			--force \
			${extra} \
			${configs} \
			--outdir "${prefix}multiqc" \
			--file-list TMP/jeter.list
		
		zip -9 -r "${prefix}multiqc.zip" "${prefix}multiqc"

##################################################################################
cat <<- EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">MultiQC</entry>
	<entry key="multiqc.version">\$( multiqc --version )</entry>
</properties>
EOF
	"""

}

