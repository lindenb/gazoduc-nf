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



/**
 * call jvarkit vcfgatkeval for one or more vcf
 *
 */
workflow JVARKIT_GATK_HARD_FILTERING_01 {
	take:
		row /* tuple [interval,vcf,vcfidx] */
	main:
		
		ch = FOR_EACH_INTERVAL(row)

		concat_ch = CONCAT_TABLES(ch.output.collect())

	emit:
		output = concat_ch.output
		pdf = concat_ch.pdf
	}

process FOR_EACH_INTERVAL {
tag "${interval} ${vcf.name}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
memory '2g'
input:
	tuple val(interval),path(vcf),path(idx)
output:
	path("*.output.table.txt"),emit:output
script:
	if(!task.ext.containsKey("percentile")) throw new IllegalArgumentException("FOR_EACH_INTERVAL percentile undefined");
	if((task.ext.percentile as double) <= 0 ) throw new IllegalArgumentException("meta.percentil <= 0");
	def percentile = task.ext.percentile as double
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP OUT


bcftools view -G --regions '${interval}' '${vcf}' |\\
	jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP vcfgatkeval --percentile ${percentile} --input-type vcf -o OUT/gatk.eval
	
mv OUT/gatk.eval.output.table.txt "${vcf.simpleName}.output.table.txt"

mv OUT TMP/
"""
}

process CONCAT_TABLES {
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
memory '2g'
input:
	path("TABLES/*")
output:
	path("OUT/gatk.eval.output.filters.txt"),emit:output
	path("OUT/gatk.eval.output.pdf"),emit:pdf
script:
	if(!task.ext.containsKey("percentile")) throw new IllegalArgumentException("FOR_EACH_INTERVAL percentile undefined");
	if((task.ext.percentile as double) <= 0 ) throw new IllegalArgumentException("meta.percentil <= 0");
	def percentile = task.ext.percentile as double

"""
hostname 1>&2
set -o pipefail
mkdir -p TMP

find TABLES  -name "*.table.txt" > TMP/jeter.list

xargs -a TMP/jeter.list -L 50 cat | jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP vcfgatkeval --percentile ${percentile} --input-type table -o TMP/gatk.eval

sed 's%output.pdf%TMP/gatk.eval.output.pdf%' TMP/gatk.eval.output.R | R --vanilla 

mv TMP OUT
"""
}
