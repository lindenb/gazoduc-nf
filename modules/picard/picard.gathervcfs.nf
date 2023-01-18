/*

Copyright (c) 2023 Pierre Lindenbaum

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

include {moduleLoad;} from '../utils/functions.nf'

process PICARD_GATHER_VCFS_01 {
tag "${list.name}"
cache "lenient"
afterScript "rm -f TMP"
cache 'lenient'
cpus 10
memory  '6g'
input:
	val(meta)
	path(list)
output:
	path("${meta.prefix?:""}gather${(meta.suffix?:".bcf").endsWith(".gz")?".vcf.gz":".bcf"}"),emit:vcf
	path("version.xml"),emit:version
	path("${meta.prefix?:""}gather${(meta.suffix?:".bcf").endsWith(".gz")?".vcf.gz.tbi":".bcf.csi"}"),emit:index
script:
	def prefix = meta.prefix?:""
	def suffix = meta.suffix?:".vcf.gz"

"""
hostname 1>&2
${moduleLoad("picard bcftools")}

mkdir TMP

java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${PICARD_JAR} GatherVcfs \
	REORDER_INPUT_BY_FIRST_VARIANT=true \
	I=${list} O=TMP/jeter.vcf.gz


if [ ! -z "${suffix.endsWith(".gz")?"":"BCF"}" ] ; then

	bcftools view --threads ${task.cpus} \
		-O ${suffix.endsWith(".gz")?"z":"b"} \
		-o "${prefix}gather${suffix}" TMP/jeter.vcf.gz

else

	mv -v TMP/jeter.vcf.gz "${prefix}gather${suffix}"

fi

bcftools index --threads ${task.cpus} ${suffix.endsWith(".gz")?"-t":""} "${prefix}gather${suffix}"

###########################################################################################"
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">Gather VCFs</entry>
	<entry key="picard">\$(java -jar \${PICARD_JAR} GatherVcfs --version 2>&1)</entry>
</properties>
EOF
"""
}
