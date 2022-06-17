
include {getKeyValue;getModules;isHg19} from '../utils/functions.nf'

process DOWNLOAD_GTF_01 {
tag "${file(reference).name}"
afterScript "rm -rf TMP"
memory "1g"
input:
	val(meta)
	val(reference)
output:
	path("${file(reference).getSimpleName()}.gtf.gz"),emit:gtf
	path("${file(reference).getSimpleName()}.gtf.gz.tbi"),emit:tbi
	path("version.xml"),emit:version
script:

	def url0 = getKeyValue(meta,"gtfurl","")
	def url = url0.isEmpty()?(isHg19(reference)?"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh37_mapping/gencode.v40lift37.annotation.gtf.gz":""):url0
	"""
	hostname 1>&2
	module load ${getModules("htslib jvarkit")}
	set -o pipefail
	mkdir TMP

	test ! -z "${url}"

	wget -O TMP/jeter0.gtf "${url}"

	if [[ `file TMP/jeter0.gtf | grep gzip` ]] ; then
		mv TMP/jeter0.gtf TMP/jeter0.gtf.gz
	fi
	

	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/bedrenamechr.jar --column 1 --convert SKIP -R "${reference}" TMP/jeter0.* > TMP/jeter.gtf
	rm TMP/jeter0.*
	mv TMP/jeter.gtf TMP/jeter0.gtf


	
	LC_ALL=C sort -S ${task.memory.kilo} -T TMP -k1,1 -k4,4n TMP/jeter0.gtf > TMP/jeter.gtf
	rm TMP/jeter0.gtf

	bgzip TMP/jeter.gtf
	tabix -p gff TMP/jeter.gtf.gz
	mv TMP/jeter.gtf.gz "${file(reference).getSimpleName()}.gtf.gz"
	mv TMP/jeter.gtf.gz.tbi "${file(reference).getSimpleName()}.gtf.gz.tbi"
	
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Download GTF for Reference</entry>
		<entry key="reference">${reference}</entry>
		<entry key="url">${url}</entry>
		<entry key="bedrenamechr">\$(java  -jar \${JVARKIT_DIST}/bedrenamechr.jar --version)</entry>
		<entry key="bgzip">\$(bgzip --version | head -n 1)</entry>
		<entry key="tabix">\$(tabix --version | head -n 1)</entry>
	</properties>
	EOF
	"""
	}
