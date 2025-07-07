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
include { SAMTOOLS_SAMPLES} from '../../samtools/samtools.samples.03.nf'
include { COLLECT_TO_FILE_01} from '../../../modules/utils/collect2file.01.nf'
include { getKeyValue;moduleLoad; assertFileExists;getVersionCmd} from '../../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.02.nf'



workflow INDEXCOV {
     take:
        fasta
	fai /* file containing the path to the bam/cram files . One per line */
	samplesheet
     main:
        ch_version = Channel.empty()

	sheet_ch  = Channel.fromPath(samplesheet).splitCsv(header:true,sep:'\t').
		branch {
			bai_only : !it.containsKey("bam") || it.get("bam").isEmpty() || it.get("bam").equals(".")
			bam_only : true
		}

	
	/** BAI only part */
	bai_only_ch = FROM_BAI_ONLY(fasta,fai, sheet_ch.bai_only.map{[it.sample,file(it.bai)]})
	bams_only_ch = REBUILD_BAI(fasta,fai, sheet_ch.bam_only.map{[it.sample,file(it.bam)]})

	

	def collate_size = ((params.batch_size?:10000) as int)

	bams2_ch = bams_only_ch.output.
		mix(bai_only_ch.output).
		toSortedList({a,b->a[0].name.compareTo(b[0].name)}).
		flatten().
		collate(collate_size).
		map{T->["batch."+T.collect{V->V.toString()}.join(" ").md5().substring(0,7) , T]}
		
	executable_ch = DOWNLOAD_GOLEFT()
	ch_version = ch_version.mix( executable_ch.version)
    	indexcov_ch = RUN_GOLEFT_INDEXCOV(fasta,fai,executable_ch.executable, bams2_ch)
	ch_version = ch_version.mix( indexcov_ch.version)

	mergebed_ch= MERGE_BEDS(indexcov_ch.output.map{T->T[0]}.collect())

	ch_version = MERGE_VERSION("Indexcov",ch_version.collect())		


	emit:
		files = indexcov_ch.files
		zip = indexcov_ch.zip
		version = ch_version
	}

process DOWNLOAD_GOLEFT {
	errorStrategy "retry"
	maxRetries 3
	output:
		path("goleft"),emit:executable
		path("version.xml"),emit:version
	script:
		def version = task.ext.version?:"v0.2.6"
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
    tag "${sample}"
    afterScript "rm -rf TMP"
    memory "5g"
    input:	
	path(fasta)
	path(fai)
	tuple val(sample),path(bai)
    output:
	path("OUT/${sample}*"),emit: output

script:
"""
	hostname 1>&2
	${moduleLoad("samtools")}

	mkdir -p TMP
	samtools dict ${fasta} > TMP/tmp.sam
	echo "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}" >> TMP/tmp.sam
	samtools sort --reference "${fasta}" -O BAM -o TMP/${sample}.bam -T TMP/tmp TMP/tmp.sam
	rm TMP/tmp.sam
	cp -v "${bai}" TMP/${sample}.bam.bai

	#just test
	samtools idxstat TMP/${sample}.bam 1>&2
	
	mv -v TMP OUT
"""
}


/**

 see https://twitter.com/yokofakun/status/1527419449669734426

 **/
process REBUILD_BAI {
    tag "${sample} ${bam.name} mapq=${params.mapq}"
    label "process_single"
    afterScript "rm -rf TMP"
    input:
        path(fasta)
        path(fai)
        tuple val(sample),path(bam)
    output:
	path("OUT/${sample}*"), emit: output
        path "version.xml",           emit: version
    script:
	def mapq = params.mapq?:1
    """
	hostname 1>&2
	${moduleLoad("samtools")}
	mkdir -p TMP

        # write header only
	samtools view --header-only -O BAM  \
		--threads ${task.cpus} \
		-o "TMP/${sample}.bam" \
		--reference "${fasta}" \
		"${bam}"

	# hack about creating bam index. see main samtools manual
	samtools view -F 3844 -q "${mapq}" --uncompressed \
		--threads ${task.cpus} \
		-o "/dev/null##idx##TMP/${sample}.bam.bai" \
		--write-index  -O BAM  \
		--reference "${fasta}" \
		"${bam}"

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
    tag "${name}"
    afterScript "rm -rf TMP"
    label "process_high"
    input:
	path(fasta)
	path(fai)
	val executable
	tuple val(name),path("BAMS/*")
    output:
    	path "*indexcov.zip"       , emit: zip
    	path "*files.list"          , emit: files
	tuple path("*.indexcov.bed.gz"), path("*.indexcov.bed.gz.tbi"),emit:output
    	path "version.xml"        , emit: version

    script:
    	def prefix = params.prefix.replaceAll("[\\.]+\$","")+"."+name
    	def prefix2 = prefix +"."
    """
	hostname 1>&2
	${moduleLoad("htslib")}
	set -o pipefail
	mkdir -p "${prefix}" TMP

	find BAMS  -name "*.bam" -o -name "*.cram" > TMP/jeter.list

	# test not empty
	test -s TMP/jeter.list

	# test all files exist
	xargs -a TMP/jeter.list -L1 --verbose test -f

	sed 's/\\.cram/.cram.crai/' TMP/jeter.list > TMP/tmp.bams.list

	# test all files exist
	xargs -a TMP/tmp.bams.list -L1 --verbose test -f


	${executable} indexcov \
		--fai "${fasta}.fai"  \
		--excludepatt `awk -F '\t' '(!(\$1 ~ /^(chr)?[0-9XY]+\$/)) {print \$1;}' "${fasta}.fai" | paste -s -d '|' `  \
		--sex `awk -F '\t' '(\$1 ~ /^(chr)?[XY]\$/) {print \$1;}' "${fasta}.fai" | paste -s -d, ` \
		--directory "${prefix}" \
		`awk '/.crai\$/ {X=1;} END {if(X==1) printf(" --extranormalize ");}' TMP/tmp.bams.list` \
		`cat TMP/tmp.bams.list`


	#create tabix index
	tabix -p bed "${prefix}/${prefix}-indexcov.bed.gz"

	ln -s "${prefix}/${prefix}-indexcov.bed.gz" "${prefix}.indexcov.bed.gz"
	ln -s "${prefix}/${prefix}-indexcov.bed.gz.tbi" "${prefix}.indexcov.bed.gz.tbi"


	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">run go left indexcov</entry>
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

process MERGE_BEDS {
tag "${L.size()}"
afterScript "rm -rf TMP"
input:
	val(L)
output:
	tuple path("${params.prefix?:""}indexcov.merged.bed.gz"),path("${params.prefix?:""}indexcov.merged.bed.gz.tbi"),emit:output
script:
"""
hostname 1>&2
${moduleLoad("htslib")}
mkdir -p TMP

# when I test in local
rm -f TMP/jeter.bed

for F in ${L.join(" ")}
do
	if ! test -f TMP/jeter.bed
	then

		gunzip -c "\${F}" | head -n 1 > TMP/header.txt
		gunzip -c "\${F}" | tail -n +2 | cut -f 1,2,3 > TMP/signature1.txt
		gunzip -c "\${F}" | tail -n +2 > TMP/jeter.bed

	else
		gunzip -c "\${F}" | tail -n +2 | cut -f 1,2,3 > TMP/signature2.txt
		# check all files have the same intervals in the same order
		cmp TMP/signature1.txt TMP/signature2.txt

		paste TMP/header.txt <(gunzip -c "\${F}" | head -n 1 | cut -f4-)  > TMP/jeter2.txt
		mv TMP/jeter2.txt TMP/header.txt

		paste TMP/jeter.bed <(gunzip -c "\${F}" | tail -n +2 | cut -f4-) > TMP/jeter2.txt
		mv TMP/jeter2.txt TMP/jeter.bed
	fi
done

cat TMP/header.txt TMP/jeter.bed > TMP/jeter2.txt
mv TMP/jeter2.txt TMP/jeter.bed

# check there is only one number of cols
test \$(awk '{print NF}' TMP/jeter.bed | uniq | sort | uniq | wc -l) -eq 1

mv TMP/jeter.bed "TMP/${params.prefix?:""}indexcov.merged.bed"


bgzip "TMP/${params.prefix?:""}indexcov.merged.bed"

tabix -p bed "TMP/${params.prefix?:""}indexcov.merged.bed.gz"

mv TMP/*.gz ./
mv TMP/*.gz.tbi ./
"""
}
