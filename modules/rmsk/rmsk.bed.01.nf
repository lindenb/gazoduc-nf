nextflow.enable.dsl=2


include {isHg19;isHg38;hasFeature;getModules} from '../utils/functions.nf'


process RMSK_TO_BED_01 {
input:
	val(meta)
	val(reference)
output:
	tuple val("rmsk"),path("rmsk.bed.gz"),path("rmsk.header"),emit:bed
	path("rmsk.bed.gz.tbi"), emit:tbi
	path("version.xml"), emit:version
script:
	def url =	(isHg19(reference)?"https://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz":
			(isHg38(reference)?"https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz":
			""));
"""
hostname 1>&2
set -o pipefail
module load ${getModules("jvarkit tabix")}

if [ ! -z "${url}" ] ; then
wget -O - "${url}" |\
	gunzip -c |\
	cut -f6-8,12 |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP |\
		LC_ALL=C sort -t '\t' -T . -k1,1 -k2,2n > rmsk.bed

else
	touch rmsk.bed
fi

bgzip rmsk.bed
tabix -p bed rmsk.bed

########################################################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<dl id="${task.process}">
	<dt>name</dt><dd>${task.process}</dd>
	<dt>description</dt><dd>repeat masker regions</dd>
	<dt>url</dt><dd>${url}</dd>
<properties>
EOF
"""
}
