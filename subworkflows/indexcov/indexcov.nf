/*

Copyright (c) 2024 Pierre Lindenbaum

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
include { SAMTOOLS_SAMPLES} from '../samtools/samtools.samples.03.nf'
include { COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include { getKeyValue;moduleLoad; assertFileExists;getVersionCmd} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'



workflow INDEXCOV {
     take:
        meta /* meta */
        genomeId
        bams /* file containing the path to the bam/cram files . One per line */
	bai_samplesheet /** tsv file containing sample and bai for BAI only files */
     main:
        ch_version = Channel.empty()


    	if ( (params.mapq  as int)> 0 || !bai_samplesheet.name.equals("NO_FILE") ) {
 		log.info("mapq ${params.mapq} is greater than 0. rebuilding bam indexes... ${bams} ${genomeId}")
		bams_ch = SAMTOOLS_SAMPLES([:], bams)
		ch_version = ch_version.mix( bams_ch.version)

        	sample_bam_fasta_ch = bams_ch.rows.filter{T->T.genomeId.equals(genomeId)}
        
		rebuild_bai_ch = REBUILD_BAI([:], sample_bam_fasta_ch)
		ch_version = ch_version.mix( rebuild_bai_ch.version.first() ) 

		if( !bai_samplesheet.name.equals("NO_FILE") ) {
			bai_only0 = FROM_BAI_ONLY(genomeId,Channel.fromPath(bai_samplesheet).splitCsv(header:true,sep:'\t'))
			bai_only = bai_only0.bam
			}
		else
			{
			bai_only = Channel.empty()
			}



		bams2_ch = COLLECT_TO_FILE_01([:], rebuild_bai_ch.bam.mix(bai_only).collect()).output
        	}
    	else
		{
		bams2_ch = Channel.fromPath(bams)
		}


	executable_ch = DOWNLOAD_GOLEFT([:])
	ch_version = ch_version.mix( executable_ch.version)
    	indexcov_ch = RUN_GOLEFT_INDEXCOV([:], executable_ch.executable ,genomeId , bams2_ch, ch_version.collect())
	ch_version = ch_version.mix( indexcov_ch.version)

	ch_version = MERGE_VERSION("Indexcov",ch_version.collect())		


	emit:
		files = indexcov_ch.files
		zip = indexcov_ch.zip
		version = ch_version
	}

process DOWNLOAD_GOLEFT {
	tag "${params.goleft.version}"
	errorStrategy "retry"
	maxRetries 3
	input:
		val(meta)
	output:
		path("goleft"),emit:executable
		path("version.xml"),emit:version
	script:
		def version = params.goleft.version
		def url ="https://github.com/brentp/goleft/releases/download/${version}/goleft_linux64"
	"""
	wget -O goleft "${url}"
	chmod +x goleft


	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Download goleft</entry>
		<entry key="version">${version}</entry>
		<entry key="url"><a>${url}</a></entry>
	</properties>
	EOF
	"""
	}

process FROM_BAI_ONLY {
    tag "${row.sample}"
    afterScript "rm -rf TMP"
    memory "5g"
    input:	
	val(genomeId)
	tuple val(row)
    output:
	path("OUT/${row.sample}.bam"),     emit: bam
	path("OUT/${row.sample}.bam.bai"), emit: bai

script:
	def genomeId = row.genomeId
        def reference = params.genomes[genomeId].fasta
	def sample  = row.sample
	def dict = "${file(reference.toRealPath()).getParent()}/${reference.getSimpleName()}.dict"
"""
	hostname 1>&2
	${moduleLoad("samtools")}

	mkdir -p TMP
	test -s "${dict}"
	cat "${dict}" > TMP/tmp.sam
	echo "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}" >> TMP/tmp.sam
	samtools sort --reference "${reference}" -O BAM -o TMP/${sample}.bam -T TMP/tmp TMP/tmp.sam
	rm TMP/tmp.sam
	cp -v "${row.bai}" TMP/${sample}.bam.bai
	mv TMP OUT	
"""
}


/**

 see https://twitter.com/yokofakun/status/1527419449669734426

 **/
process REBUILD_BAI {
    tag "${row.sample} ${file(row.bam).name} mapq=${params.mapq}"
    afterScript "rm -rf TMP"
    memory "5g"
    cpus 3
    input:
        val meta
        val(row)
    output:
	path("OUT/${row.sample}.bam"),     emit: bam
	path("OUT/${row.sample}.bam.bai"), emit: bai
        path "version.xml",           emit: version

    script:
	def genomeId = row.genomeId
	def reference = params.genomes[genomeId].fasta
	def mapq = params.mapq?:1
    """
	hostname 1>&2
	${moduleLoad("samtools")}
	mkdir -p TMP

        # write header only
	samtools view --header-only -O BAM  \
		--threads ${task.cpus} \
		-o "TMP/${row.sample}.bam" \
		--reference "${reference}" \
		"${row.bam}"

	# hack about creating bam index. see main samtools manual
	samtools view -F 3844 -q "${mapq}" --uncompressed \
		--threads ${task.cpus} \
		-o "/dev/null##idx##TMP/${row.sample}.bam.bai" \
		--write-index  -O BAM  \
		--reference "${reference}" \
		"${row.bam}"

	mv -v TMP OUT

	##########################################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Re-Build BAI without generating a BAM</entry>
		<entry key="mapq">${mapq}</entry>
		<entry key="samtools">${getVersionCmd("samtools")}</entry>
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
	val genomeId
	path bams
	val(versions)
    output:
    	path "*indexcov.zip"       , emit: zip
    	path "*files.list"          , emit: files
    	path "version.xml"        , emit: version

    script:
	def fasta = params.genomes[genomeId].fasta
    	def prefix = params.prefix.replaceAll("[\\.]+\$","")
    	def prefix2 = prefix +"."
    """
	hostname 1>&2
	${moduleLoad("htslib")}
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
		<entry key="versions">${getVersionCmd("tabix")}</entry>
	</properties>
	EOF

	# all the generated files
	find "${prefix}" -type f > "${prefix2}files.list"

	# zip results
	cat "${prefix2}files.list" | zip -@ -9  "${prefix2}indexcov.zip"
    """
}

