include { SAMTOOLS_SAMPLES01} from '../../modules/samtools/samtools.samples.01.nf'
include { COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include { getKeyValue; getModules; assertFileExists} from '../../modules/utils/functions.nf'

workflow INDEXCOV {
     take:
        meta /* meta */
        reference /* indexed fasta reference */
        bams /* file containing the path to the bam/cram files . One per line */
     main:
	assertFileExists(reference,"reference must be defined")
	assertFileExists(bams,"path to bams must be defined")
        mapq = ((meta.mapq?:"0") as int)
        ch_version = Channel.empty()


    	if ( mapq > 0) {
 		log.info("mapq ${mapq} is greater than 0. rebuilding bam indexes... ${bams} ${reference}")
		bams_ch = SAMTOOLS_SAMPLES01([:], reference, bams)
		ch_version = ch_version.mix( bams_ch.version)

        	sample_bam_fasta_ch = bams_ch.output.splitCsv(header: false,sep:'\t',strip:true)
        
		rebuild_bai_ch = REBUILD_BAI(meta.subMap(["mapq"]), reference, sample_bam_fasta_ch)
		ch_version = ch_version.mix( rebuild_bai_ch.version.first() ) 

		bams2_ch = COLLECT_TO_FILE_01([:], rebuild_bai_ch.bam.collect()).output
        	}
    	else
		{
		bams2_ch = Channel.fromPath(bams)
		}

	executable_ch = DOWNLOAD_GOLEFT(meta.subMap(["goleft_version"]))
	ch_version = ch_version.mix( executable_ch.version)
    	indexcov_ch = RUN_GOLEFT_INDEXCOV(meta, executable_ch.executable ,reference , bams2_ch, ch_version.collect())

	emit:
		files = indexcov_ch.files
		zip = indexcov_ch.zip
		version = indexcov_ch.version
	}

process DOWNLOAD_GOLEFT {
	errorStrategy "retry"
	maxRetries 3
	input:
		val(meta)
	output:
		path("goleft"),emit:executable
		path("version.xml"),emit:version
	script:
		def version = getKeyValue(meta,"goleft_version","v0.2.4")
		def url ="https://github.com/brentp/goleft/releases/download/${version}/goleft_linux64"
	"""
	wget -O goleft "${url}"
	chmod +x goleft


	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Download goleft</entry>
		<entry key="version">${version}</entry>
		<entry key="url">${url}</entry>
	</properties>
	EOF
	"""
	}


/**

 see https://twitter.com/yokofakun/status/1527419449669734426

 **/
process REBUILD_BAI {
    tag "${sample} ${file(bam).name} mapq=${meta.mapq}"
    afterScript "rm -rf TMP"
    label "process_low"
    input:
        val meta
	val reference
        tuple val(sample), val(bam)
    output:
	path("OUT/${sample}.bam"),     emit: bam
	path("OUT/${sample}.bam.bai"), emit: bai
        path "version.xml",           emit: version

    script:
	def mapq = getKeyValue(meta,"mapq","1")
    """
	mkdir TMP

        # write header only
	samtools view --header-only -O BAM  \
		--threads ${task.cpus} \
		-o "TMP/${sample}.bam" \
		--reference "${reference}" \
		"${bam}"

	# hack about creating bam index. see main samtools manual
	samtools view -F 3844 -q "${mapq}" --uncompressed \
		--threads ${task.cpus} \
		-o "/dev/null##idx##TMP/${sample}.bam.bai" \
		--write-index  -O BAM  \
		--reference "${reference}" \
		"${bam}"

	mv TMP OUT

	##########################################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Re-Build BAI without generating a BAM</entry>
		<entry key="mapq">${mapq}</entry>
		<entry key="samtools">\$(samtools  --version | head -n 1| cut -d ' ' -f2)</entry>
	</properties>
	EOF
    """
   }


process RUN_GOLEFT_INDEXCOV {
    tag "${bams.name}"
    label "process_high"
    input:
	val meta
	val executable
	val fasta
	path bams
	val(versions)
    output:
    	path "*indexcov.zip"       , emit: zip
    	path "*files.list"          , emit: files
    	path "version.xml"        , emit: version

    script:
    	def prefix = getKeyValue(meta,"prefix","indexcov").replaceAll("[\\.]+\$","")
    	def prefix2 = prefix +"."
    """
	hostname 1>&2
	module load ${getModules("htslib")}
	set -o pipefail
	mkdir -p "${prefix}"

	# test not empty
	test -s "${bams}"

	# test all files exist
	xargs -a "${bams}" -L1 --verbose test -f

	sed 's/\\.cram/.cram.crai/' "${bams}" > tmp.bams.list

	# test all files exist
	xargs -a tmp.bams.list -L1 --verbose test -f

	

	${executable} indexcov \
		--fai "${fasta}.fai"  \
		--excludepatt `awk -F '\t' '(!(\$1 ~ /^(chr)?[0-9XY]+\$/)) {print \$1;}' "${fasta}.fai" | paste -s -d '|' `  \
		--sex `awk -F '\t' '(\$1 ~ /^(chr)?[XY]\$/) {print \$1;}' "${fasta}.fai" | paste -s -d, ` \
		--directory "${prefix}" \
		`awk '/.crai\$/ {X=1;} END {if(X==1) printf(" --extranormalize ");}' tmp.bams.list` \
		`cat tmp.bams.list`

	rm tmp.bams.list

	#create tabix index
	tabix -p bed "${prefix}/${prefix}-indexcov.bed.gz"




	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">run go left indexcov</entry>
		<entry key="bams">${bams}</entry>
		<entry key="count-bams">\$(wc -l < {bams})</entry>
		<entry key="goleft">\$(${executable} -h |  head -n 1 | sed 's/.*: //')</entry>
		<entry key="steps">
	EOF

	for X in ${versions.join(" ")}
	do
		xmllint --format "\${X}" | tail -n+2 >> version.xml
	done

	cat <<- EOF >> version.xml
		</entry>
	</properties>
	EOF
	xmllint --format version.xml > "${prefix}/${prefix}-version.xml"

	# all the generated files
	find "${prefix}" -type f > "${prefix2}files.list"

	# zip results
	cat "${prefix2}files.list" | zip -@ -9  "${prefix2}indexcov.zip"

    """
}

