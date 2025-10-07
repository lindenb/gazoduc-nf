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


process BED_CLUSTER {
	label "process_single"
	tag "${meta.id?:""} "
	afterScript "rm -rf TMP"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	when:
	    task.ext.when == null || task.ext.when
	input:
		tuple val(meta1),path(fasta)
		tuple val(meta2),path(fai)
		tuple val(meta3),path(dict)
		tuple val(meta ),path(bed)
	output:
		tuple val(meta), path("BEDS/*.bed",arity:"0..*"),emit:bed
		path("versions.yml"),emit:versions
	script:
		def args = task.ext.args?:""
		if(args.trim().isEmpty()) throw new IllegalArgumentException("args for bedcluter must be defined in ${task.process}")
	"""
	hostname 1>&2
	set -o pipefail
	mkdir -p TMP/TMP2 BEDS
	jvarkit  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData bedcluster \\
		-R "${fasta}" \\
		${args} \\
		-o TMP/TMP2 "${bed}"

	#force uniq names for later
	MD5=\$(pwd | sha1sum  | cut -c 1-10)
	find TMP/TMP2 -type f -name "*.bed" |\
	while read F
	do
		mv "\${F}" "BEDS/\${MD5}.\$(basename "\$F")"
	done

cat << EOF > versions.yml
${task.process}:
	jvarkit: "\$(jvarkit --version)"
EOF
	"""
	
stub:
"""
mkdir BEDS
split --additional-suffix=.bed --lines=1 "${bed}" BEDS/cluster
touch versions.yml
"""	
}
