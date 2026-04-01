/*

Copyright (c) 2026 Pierre Lindenbaum

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
include {parseBoolean;isBlank} from '../utils/functions.nf'

process MOSDEPTH_RUN_01 {
	tag "${row.sample} ${file(row.bam).name}"
	label "process_single"
	afterScript "rm -rf TMP"
	input:
		path(executable)
		path(reference)
		tuple val(sample),path(bam),path(bai),path(bed)
	output:
		tuple val(sample),path("${sample}.*"),emit:output
	script:
		def fasta = reference.find{it.name.endsWith("a")}.first()
		def mapq = task.ext.mapq
		def suffix = sample
		def extra = task.ext.args
	"""
	hostname 1>&2
	mkdir -p TMP

        ${executable} ${bed.name.equals("NO_FILE")?"":"--by \"${bed}\""} \
 		-t ${task.cpus} --fasta "${fasta}" --mapq ${mapq} ${extra} \
		'TMP/${suffix}' "${bam}"

	mv -v TMP/${suffix}* ./
	"""
	}
