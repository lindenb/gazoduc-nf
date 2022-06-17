include {getKeyValue;getModules} from '../../modules/utils/functions.nf'

process GTF_TO_UORF_01 {
tag "${file(gtf).name}"
afterScript "rm -rf TMP"
memory "5g"
input:
	val(meta)
	val(reference)
	val(gtf)
output:
	path("uorf.gtf.gz"),emit:gtf
	path("uorf.gtf.gz.tbi"),emit:tbi
	path("version.xml"),emit:version
script:
	def strength = getKeyValue(meta,"strength","Strong")
"""
hostname 1>&2
set -o pipefail
module load ${getModules("jvarkit htslib")}

mkdir TMP

(gunzip -c "${gtf}" || cat "${gtf}" ) |\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/gtfupstreamorf.jar -R "${reference}" --strength "${strength}" |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k4,4n |\
	bgzip > uorf.gtf.gz

tabix -p gff uorf.gtf.gz

#################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}<entry>
	<entry key="name">extract uORF <entry>
	<entry key="input">${gtf}<entry>
	<entry key="strength">${strength}<entry>
	<entry key="reference">${reference}<entry>
        <entry key="gtfupstreamorf">\$(java  -jar \${JVARKIT_DIST}/gtfupstreamorf.jar --version)</entry>
        <entry key="bgzip">\$(bgzip --version | head -n 1)</entry>
        <entry key="tabix">\$(tabix --version | head -n 1)</entry>
</properties>
EOF
"""
}
