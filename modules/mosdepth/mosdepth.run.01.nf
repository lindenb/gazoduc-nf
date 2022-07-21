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
include {parseBoolean;isBlank} from '../utils/functions.nf'

process MOSDEPTH_RUN_01 {
	tag "${row.sample} ${row.bam}"
	cpus 1
	input:
		val(meta)
		val(executable)
		val(row)
	output:
		tuple val(row),path("${row.sample}${row.suffix?:""}.mosdepth.global.dist.txt"),emit:globaldist
		tuple val(row),path("${row.sample}${row.suffix?:""}.mosdepth.summary.txt"),emit:summary
		tuple val(row),path("${row.sample}${row.suffix?:""}.per-base.bed.gz"),optional:true,emit:perbase
		path("version.xml"),emit:version
	script:
		def bed = row.bed?:""
		def median = parseBoolean(row.use_median?:"true")?"--use-median":""
		def per_base = parseBoolean(row.per_base?:"false")?"--no-per-base":""
		def mapq =row.mapq?:"0"
		def suffix = row.suffix?:""
	"""
	hostname 1>&2
	

        ${executable} ${isBlank(bed)?"":"--by \""+bed+"\""} ${median} ${per_base} \
 		-t ${task.cpus} --fasta "${row.reference}" --mapq ${mapq} \
		${row.sample}${suffix} ${row.bam}



##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">apply mosdepth to bam</entry>
	<entry key="bam">${row.bam}</entry>
        <entry key="mosdepth.version">\$(${executable} --version)</entry>
        <entry key="median">${median}</entry>
        <entry key="mapq">${mapq}</entry>
</properties>
EOF
	"""
	}
