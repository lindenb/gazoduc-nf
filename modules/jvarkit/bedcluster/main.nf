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


process BED_CLUSTER {
	label "process_single"
	tag "${meta.id?:""} "
	afterScript "rm -rf TMP"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	input:
		tuple val(meta1),path(dict)
		tuple val(meta ),path(bed)
	output:
		tuple val(meta), path("*.bed",arity:"0..*"),emit:bed
		path("versions.yml"),emit:versions
	script:
		def args = task.ext.args?:""
		if(args.trim().isEmpty()) throw new IllegalArgumentException("args for bedcluter must be defined in ${task.process}")
		def prefix = task.ext.prefix?:"${meta.id}"
		def jvm = task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData"
		def jvarkit  = "java ${jvm} -jar \${HOME}/jvarkit.jar" //todo , use conda version
	"""
	mkdir -p TMP/TMP2
	${jvarkit} bedcluster \\
		-R "${dict}" \\
		--prefix "${prefix}" \\
		${args} \\
		-o TMP/TMP2 "${bed}"

	#  move files if they exists
	mv -v TMP/TMP2/*.bed ./ || true

cat << EOF > versions.yml
${task.process}:
	jvarkit: "\$(jvarkit --version)"
EOF
"""
	
stub:
	def prefix = task.ext.prefix?:"${meta.id}"
"""
mkdir BEDS
split --additional-suffix=.bed --lines=1 "${bed}" BEDS/${prefix}.cluster
mv BEDS/*.bed ./ || true
touch versions.yml
"""	
}
