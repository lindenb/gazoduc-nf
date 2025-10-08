/*

Copyright (c) 2025 Pierre Lindenbaum

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

process EXTRACT_BAM {
	tag "${meta.id?:bam.name}"
	label "process_single"
	conda "${moduleDir}/../../../conda/somalier.yml"

	input:
		tuple val(meta1),path(sites),path(sites_idx)
		tuple val(meta),path(bam),path(bai),path(fasta),path(fai)
	output:
		tuple val(meta),path("extracted/*.somalier"),emit:output
		path("versions.yml"),emit:versions
	script:
		def prefix = meta.id?:"${bam.name}"
	"""
	hostname 1>&2
	mkdir -p extracted
	somalier extract -d extracted --sites "${sites}" -f "${fasta}" "${bam}"
	
	test -s extracted/*.somalier

	# change name if needed
	if ! test -f extracted/${prefix}.somalier
	then
		mv  -v extracted/*.somalier "extracted/${prefix}.somalier"
	fi


cat << EOF > versions.yml
${task.process}:
    somalier: todo
EOF
	"""
stub:
"""
mkdir extracted
touch versions.yml "extracted/${meta.id}.somalier" 
"""
}
