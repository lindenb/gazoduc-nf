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
nextflow.enable.dsl=2

/** path to indexed fasta reference */
params.reference = ""
params.vcf = "NO_FILE"
params.bed = ""
params.samplesheet = "NO_FILE"
params.publishDir = ""
params.prefix = ""

include {VG_CONSTRUCT_01} from '../../subworkflows/vg/vg.construct.01.nf'
include {VG_DOWNLOAD_01} from '../../modules/vg/vg.download.01.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {runOnComplete} from '../../modules/utils/functions.nf'


workflow {
	ch1 = VG_MAP_FASTQS(params,params.reference, Channel.fromPath(params.samplesheet), file(params.bed), file(params.vcf))
	html = VERSION_TO_HTML(params,ch1.version)
	//PUBLISH(ch1.version,html.html,ch1.vcf,ch1.index,ch1.pdf)
	}
runOnComplete(workflow);

workflow VG_MAP_FASTQS {
	take:
		meta
		reference
		samplesheet
		bed
		vcf
	main:
		version_ch = Channel.empty()

		ex_ch = VG_DOWNLOAD_01([:])
		version_ch = version_ch.mix(ex_ch.version)

		vg_ch = VG_CONSTRUCT_01([:], reference , ex_ch.executable, vcf, bed)
		version_ch = version_ch.mix(vg_ch.version)

		each_sample = samplesheet.splitCsv(sep:',',header:true).
				combine(vg_ch.vg).
				map{T->[[:],T[0].plus(T[1])]}


		gam_ch = VG_MAP(ex_ch.executable, each_sample)

		to_call_ch = gam_ch.output.map{T->[T[0],T[1].plus(gam:T[2])]}
	
		VG_CALL(ex_ch.executable,to_call_ch)
	
	emit:
		version = version_ch	
	}

process VG_MAP {
	tag "${row.sample}"
	input:
		val(executable)
		tuple val(meta),val(row)
	output:
		tuple val(meta),val(row),path("${row.sample}.gam"),emit:output
	script:
	"""
	hostname 1>&2
	mkdir -p TMP
	echo '${row}'

	export TMPDIR=\${PWD}/TMP

	${executable} map --fastq "${row.R1}" --fastq "${row.R2}" -x "${row.xg}" -g "${row.gcsa}" > TMP/jeter.gam

	mv -v TMP/jeter.gam "${row.sample}.gam"

	"""
	}

process VG_CALL {
	tag "${row.sample}"
	afterScript "rm -rf TMP"
	input:
		val(executable)
		tuple val(meta),val(row)
	output:
		path("${row.sample}.bcf"),emit:vcf
		path("${row.sample}.bcf.csi"),emit:index
	script:
	"""
	hostname 1>&2
	mkdir -p TMP
	module load bcftools/0.0.0
	set -o pipefail
	set -x
	export TMPDIR=\${PWD}/TMP

	${executable} augment "${row.vg}"  "${row.gam}"  -A TMP/augment.gam > augment.vg
	${executable} pack -x  "${row.xg}" -g augment.vg -Q 5 -s 5 -o TMP/jeter.pack

	${executable} call "${row.xg}" --snarls "${row.snarls}" --sample "${row.sample}" -k TMP/jeter.pack -g "${row.gcsa}" | bcftools view -O b -o TMP/${row.sample}.bcf

	bcftools index  TMP/${row.sample}.bcf
	mv  TMP/${row.sample}.bcf ./
	mv  TMP/${row.sample}.bcf.csi ./
	"""
	}




