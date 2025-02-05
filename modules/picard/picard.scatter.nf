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




process SCATTER_INTERVALS_BY_NS {
label "process_quick"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	path(reference) // fasta,fai,dict
output:
	path("*.interval_list"), emit:output
script:
	def fasta = reference.find{it.name.endsWith("a")}.first()
	if(!task.ext.containsKey("type")) throw new IllegalArgumentException("SCATTER_INTERVALS_BY_NS type missing")
	if(!task.ext.containsKey("max_to_merge")) throw new IllegalArgumentException("SCATTER_INTERVALS_BY_NS type missing")
	def type = task.ext.type
	def max_to_merge = task.ext.max_to_merge
	"""
	hostname 1>&2
	mkdir -p TMP
	gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" ScatterIntervalsByNs \\
	    -R "${fasta}" \\
	    --MAX_TO_MERGE "${max_to_merge}" \\
	    -O "${fasta.getSimpleName()}.${type}.${max_to_merge}.interval_list" \\
	    -OUTPUT_TYPE "${type}"
	"""
	}

