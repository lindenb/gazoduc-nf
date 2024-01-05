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
include {moduleLoad;getVersionCmd} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'




/**
 * call jvarkit vcfgatkeval for one or more vcf
 *
 */
workflow JVARKIT_GATK_HARD_FILTERING_01 {
	take:
		meta /* contains: percentile */
		row /* contains: vcf  interval */
	main:
		if(!meta.containsKey("percentile")) throw new IllegalArgumentException("percentile undefined");
		if((meta.percentile as double) <= 0 ) throw new IllegalArgumentException("meta.percentil <= 0");

		version_ch = Channel.empty()
		
		rch = FOR_EACH_INTERVAL(meta, row)
		version_ch = version_ch.mix(rch.version)

		concat_ch = CONCAT_TABLES(meta,rch.output.collect())
		version_ch = version_ch.mix(concat_ch.version)

		version_ch = MERGE_VERSION("gatk eval hard filtering",version_ch.collect())
	emit:
		output = concat_ch.output
		pdf = concat_ch.pdf
		version = version_ch
	}

process FOR_EACH_INTERVAL {
tag "${row.interval} ${file(row.vcf).name}"
afterScript "rm -rf TMP"
memory '2g'
input:
	val(meta)
	val(row)
output:
	path("OUT/gatk.eval.output.table.txt"),emit:output
	path("version.xml"),emit:version
script:
	if(!meta.containsKey("percentile")) throw new IllegalArgumentException("percentile undefined");
	if(!row.containsKey("vcf")) throw new IllegalArgumentException("vcf undefined");
	if(!row.containsKey("interval")) throw new IllegalArgumentException("interval undefined");
	if((meta.percentile as double) <= 0 ) throw new IllegalArgumentException("meta.percentil <= 0");
"""
hostname 1>&2
${moduleLoad("jvarkit bcftools")}
set -o pipefail
mkdir -p TMP OUT


bcftools view -G --regions '${row.interval}' '${row.vcf}' |\
	java -jar -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP \${JVARKIT_DIST}/jvarkit.jar vcfgatkeval --percentile ${meta.percentile} --input-type vcf -o OUT/gatk.eval
	
	
###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="versions">${getVersionCmd("bcftools jvarkit/jvarkit")}</entry>
</properties>
EOF
"""
}

process CONCAT_TABLES {
tag "N=${L.size()}"
afterScript "rm -rf TMP"
memory '2g'
input:
	val(meta)
	val(L)
output:
	path("OUT/gatk.eval.output.filters.txt"),emit:output
	path("OUT/gatk.eval.output.pdf"),emit:pdf
	path("version.xml"),emit:version	
script:
	if(!meta.containsKey("percentile")) throw new IllegalArgumentException("percentile undefined");
	if((meta.percentile as double) <= 0 ) throw new IllegalArgumentException("meta.percentil <= 0");

"""
hostname 1>&2
${moduleLoad("jvarkit R")}
set -o pipefail
mkdir -p TMP

cat << EOF > TMP/jeter.list
${L.join("\n")}
EOF

xargs -a TMP/jeter.list -L 50 cat | java -jar -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP \${JVARKIT_DIST}/jvarkit.jar vcfgatkeval --percentile ${meta.percentile} --input-type table -o TMP/gatk.eval

sed 's%output.pdf%TMP/gatk.eval.output.pdf%' TMP/gatk.eval.output.R | R --vanilla 

mv TMP OUT
###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="versions">${getVersionCmd("bcftools jvarkit/jvarkit")}</entry>
</properties>
EOF
"""
}
