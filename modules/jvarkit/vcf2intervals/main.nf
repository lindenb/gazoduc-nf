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

process JVARKIT_VCF_TO_INTERVALS {
	label "process_single"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	afterScript "rm -rf TMP"
	tag "${interval} ${vcf.name}"
	input:
		tuple val(meta),path(vcf),path(idx),path(optional_bed)
	output:
		tuple val(meta),path("*.bed"),optional:true,emit:bed /* chrom/start/end/vcf/idx */
		path("versions.yml"),emit:versions
	script:
		if(task.ext.distance==null) throw new IllegalStateException("task.ext.distance missing for ${task.process}");
		if(task.ext.min_distance==null) throw new IllegalStateException("task.ext.min_distance missing for ${task.process}");
		def distance = task.ext.distance?:"10mb"
		def min_distance = task.ext.min_distance?:"1kb"
		def cpu2 = java.lang.Math.max(1,task.cpus-2)
	"""
	hostname 1>&2
	set -o pipefail

	mkdir -p TMP

	bcftools view --threads ${cpu2}  ${optional_bed?"--regions-file \"${optional_bed.name}\""} "${vcf}" |\\
	jvarkit  -Djava.io.tmpdir=TMP -Xmx${task.memory.giga}G vcf2intervals  \\
		--bed \\
		--distance "${distance}" \\
		--min-distance "${min_distance}" |\\
		awk  '{printf("%s\t%s\t%s\t${vcf.toRealPath()}\t${idx.toRealPath()}\\n",\$1,\$2,\$3);}' > TMP/intervals.bed

	MD5=`echo '${interval} ${bed} ${vcf} ${idx}' | sha1sum | cut -d ' ' -f1`	
	mv TMP/intervals.bed "\${MD5}.vcf2interval.bed"

cat << EOF > versions.yml
${task.process}:
	jvarkit: TODO
EOF
	"""
	}
