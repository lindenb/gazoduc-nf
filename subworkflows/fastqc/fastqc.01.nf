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
include {moduleLoad;escapeXml} from '../../modules/utils/functions.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {APPLY_FASTQC_01} from '../../modules/fastqc/fastqc.01.nf'
include {MULTIQC_01} from '../../modules/multiqc/multiqc.01.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {SIMPLE_ZIP_01} from '../../modules/utils/zip.simple.01.nf'

workflow FASTQC_01 {
	take:
		meta
		fastqs
	main:
		version_ch = Channel.empty()

		cont_ch = WGET_CONTAMINANTS(meta)
		version_ch = version_ch.mix(cont_ch.version)

		adapt_ch = WGET_ADAPTERS(meta)
		version_ch = version_ch.mix(adapt_ch.version)

		each_fastq_ch = fastqs.splitText().map{it->it.trim()}.
			filter{L->!L.startsWith("#")}.
			filter{L->!L.isEmpty()}.
			map{T->[
				fastq:T
				]}.
			combine(cont_ch.output).
			map{T->T[0].plus("contaminants":T[1])}.
			combine(adapt_ch.output).
			map{T->T[0].plus("adapters":T[1])}

		to_zip = Channel.empty()

		qc_ch = APPLY_FASTQC_01(meta,each_fastq_ch)
		version_ch = version_ch.mix(qc_ch.version)
		to_zip = to_zip.mix(qc_ch.output.map{T->T[1]})

		file_list_ch = COLLECT_TO_FILE_01([:],qc_ch.output.map{T->T[1]}.collect())
		version_ch = version_ch.mix(file_list_ch.version)
		

		multiqc_ch = MULTIQC_01(meta,file_list_ch.output.map{T->["files":T,"prefix":(meta.prefix?:"")]})
		version_ch = version_ch.mix(multiqc_ch.version)
		to_zip = to_zip.mix(multiqc_ch.zip)


		version_ch = MERGE_VERSION(meta , "FASTQC", "FASTC", version_ch.collect())
		to_zip = to_zip.mix(version_ch)

		html = VERSION_TO_HTML(meta, version_ch)
		to_zip = to_zip.mix(html.html)

		zip_ch = SIMPLE_ZIP_01(meta,to_zip.collect())

	emit:
		zip = zip_ch.zip
		version = version_ch
		multiqc = multiqc_ch.zip
	}

process WGET_CONTAMINANTS {
	executor "local"
	input:
		val(meta)
	output:
		path("contaminant_list.txt"), emit:output
		path("version.xml"),emit:version
	script:
		def url = "https://raw.githubusercontent.com/s-andrews/FastQC/master/Configuration/contaminant_list.txt"
	"""
	hostname 1>&2
	wget -O contaminant_list.txt "${url}"


##################################################################################
cat <<- EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">Download contamination file</entry>
	<entry key="wget.version"><url>${url}</url></entry>
	<entry key="wget.version">\$(wget --version | head -n1)</entry>
</properties>
EOF
	"""
	}

process WGET_ADAPTERS {
	executor "local"
	input:
		val(meta)
	output:
		path("adapter_list.txt"), emit:output
		path("version.xml"),emit:version
	script:
		def url = "https://raw.githubusercontent.com/s-andrews/FastQC/master/Configuration/adapter_list.txt"
	"""
	hostname 1>&2
	wget -O adapter_list.txt "${url}"


##################################################################################
cat <<- EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">Download adapters file</entry>
	<entry key="wget.version"><url>${url}</url></entry>
	<entry key="wget.version">\$(wget --version | head -n1)</entry>
</properties>
EOF
	"""
	}
