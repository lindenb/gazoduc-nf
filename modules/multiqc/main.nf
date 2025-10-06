

process MULTIQC {
label "process_quick"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/multiqc.yml"
input:
	tuple val(meta),path(multiqc_files, stageAs: "?/*")
output:
	tuple val(meta),path("*.zip"),emit:zip
    path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:""
    def comment = task.ext.comment ?:""
    def title = task.ext.comment ?:""
"""
mkdir -p TMP
export TMPDIR=\${PWD}/TMP

mkdir -p "${prefix}multiqc"

export LC_ALL=en_US.utf8

multiqc --filename  "${prefix}multiqc_report.html" \\
    --no-ansi \\
	--title "${title}" \\
	--comment "${comment}"  \\
	--force \\
	--outdir "${prefix}multiqc"  \\
	.
		
zip -9 -r "${prefix}multiqc.zip" "${prefix}multiqc"


cat << EOF > versions.yml
${task.process}:
    multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
EOF
"""

stub:
"""
touch multiqc.zip
touch versions.yml
"""
}