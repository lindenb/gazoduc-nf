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


process COVERAGE_GRID {
	label "process_single"
	tag "${meta.id?:""}"
	afterScript "rm -rf TMP"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	input:
		tuple val(meta1),path(fasta)
		tuple val(meta2),path(fai)
		tuple val(meta3),path(dict)
		tuple val(meta4),path(gtf),path(gtf_tbi)
		tuple val(meta5),path(known_vcf),path(known_vcf_tbi)
		tuple val(meta6),path(samplesheet)
		tuple val(meta7),path("BAMS/*")
		tuple val(meta ),val(contig),val(start),val(end)
	output:
		tuple val(meta), path("*.ps"),emit:postscript
		path("versions.yml"),emit:versions
	script:
		def mapq = task.ext.mapq?:"-1"
		def args1 = task.ext.args1?:""
		def title = task.ext.title?:"${contig}:${start}-${end}"
		def extend = task.ext.extend?:"3"
		def dimension = task.ext.dimension?:"1618x1000"
		def prefix = task.ext.prefix?:"${meta.id}.covgrid"
        def jvm =  task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData"
        def jvarkit = task.ext.jvarkit?:"java ${jvm} -jar \${HOME}/jvarkit.jar"// TODO update when conda released
	"""
	mkdir -p TMP
	if ${samplesheet?true:false}
	then
		cp ${samplesheet} TMP/bams.tsv
	else
	find BAMS/ \\( -name "*.bam" -o -name "*.cram" \\) |\\
		samtools samples |\\
		sort -T TMP -t '\t' -k1,1 |\\
		awk -F '\t' 'BEGIN{printf("sample\tbam\\n");} {print;}' > TMP/bams.tsv
	fi
	test -s TMP/bams.tsv

	${jvarkit} coveragegrid \\
		-R ${fasta} \\
		--region "${contig}:${start}-${end}" \\
		--dimension "${dimension}" \
		${(mapq as int) >= 0?"--mapq ${mapq}":""} \\
		--threads ${task.cpus} \\
		--extend ${extend} \\
		--title "${title}" \\
		--format PS \\
		${gtf?"--gtf \"${gtf}\"":""} \\
		${known_vcf?"--known \"${known_vcf}\"":""} \\
		${args1} \\
		-o "TMP/jeter.ps" \\
		${samplesheet}

mv TMP/jeter.ps ${prefix}.ps

cat << EOF > versions.yml
${task.process}:
	jvarkit: "\$(jvarkit --version)"
EOF
"""
	
stub:
def prefix = task.ext.prefix?:"${meta.id}.covgrid"
"""
touch ${prefix}.ps
touch versions.yml
"""	
}
