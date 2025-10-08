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
process JVARKIT_BAM_WITHOUT_BAI {
	label "process_single"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	afterScript "rm -rf TMP"
	tag "${meta.id?:""} ${url} ${meta.contig}:${meta.start}-${meta.end}"
	input:
		tuple val(meta1),path(dict)
		tuple val(meta ),val(url)
	output:
		tuple val(meta ),path("*.bam"),path("*.bai"),emit:bam
		path("versions.yml"),emit:versions
	script:
		def args1 = task.ext.args1?:""
		if(!(meta.contig)) throw new IllegalArgumentException("meta.contig missing");
		if(!(meta.start)) throw new IllegalArgumentException("meta.start missing");
		if(!(meta.end)) throw new IllegalArgumentException("meta.end missing");
		def prefix =  task.ext.prefix?:"${meta.id?:""}.${meta.contig}_${meta.start}_${meta.end}"
	"""
	hostname 1>&2
	mkdir -p TMP

	env | grep -i proxy 1>&2

	java  -XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar ${HOME}/jvarkit.jar bamwithoutbai \\
		-r "${meta.contig}:${meta.start}-${meta.end}"  ${args1}  '${url}' |\\
		samtools addreplacerg -r '@RG\\tID:${meta.id}\\tSM:${meta.id}' -o TMP/jeter.bam -O BAM -
	
	if ${dict?true:false}
	then

		java  -XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar ${HOME}/jvarkit.jar bamrenamechr \\
			--dict "${dict}" \\
			--samoutputformat BAM \\
			-o TMP/jeter2.bam \\
			TMP/jeter.bam

		mv -v TMP/jeter2.bam TMP/jeter.bam

	fi

	samtools index -@ ${task.cpus} TMP/jeter.bam

	mv -v TMP/jeter.bam ${prefix}.bam
	mv -v TMP/jeter.bam.bai ${prefix}.bam.bai

cat << EOF > versions.yml
${task.process}:
	jvarkit: "\$(jvarkit --version)"
EOF
"""
stub:
def prefix =  task.ext.prefix?:"${meta.id?:""}.${meta.contig}_${meta.start}_${meta.end}"
"""
touch versions.yml ${prefix}.bam ${prefix}.bam.bai
"""
}
