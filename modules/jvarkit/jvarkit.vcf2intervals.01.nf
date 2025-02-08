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



process JVARKIT_VCF_TO_INTERVALS_01 {
	label "process_short"
	conda "${moduleDir}/../../conda/bioinfo.01.yml"
	afterScript "rm -rf TMP"
	tag "${interval} ${vcf.name}"
	input:
		tuple val(interval),path(bed),path(vcf),path(idx)
	output:
		path("*.bed"),emit:bed /* chrom/start/end/vcf */
	script:
		def distance = task.ext.distance?:"10mb"
		def min_distance = task.ext.min_distance?:"1kb"
	"""
	hostname 1>&2
	set -o pipefail

	mkdir -p TMP

	bcftools view -G "${vcf}" ${bed.name.endsWith(".bed")?"\"${bed}\"":"\"${interval}\""} |\\
	jvarkit  -Djava.io.tmpdir=TMP -Xmx${task.memory.giga}G vcf2intervals  \\
		--bed \\
		--distance "${distance}" \\
		--min-distance "${min_distance}" |\\
		awk  '{printf("%s\t%s\t%s\t${vcf.toRealPath()}\t${idx.toRealPath()}\\n",\$1,\$2,\$3);}' > TMP/intervals.bed

	MD5=`echo '${interval} ${bed} ${vcf} ${idx}' | sha1sum | cut -d ' ' -f1`	
	mv TMP/intervals.bed "\${MD5}.vcf2interval.bed"
	"""
	}
