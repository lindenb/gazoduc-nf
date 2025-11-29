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
include { isBlank      } from '../../utils/functions.nf'
include { verify       } from '../../utils/functions.nf'


process GTF_TO_EXOME {
tag "${meta.id}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(fai) //for slop
	tuple val(meta),path(gtf)
output:
	tuple val(meta),path("*.bed"),emit:bed
	path("versions.yml"),emit:versions
script:
	def what = task.ext.what?:""
    verify(!isBlank(what),"${task.process} : task.ext.what must be defined (e.g: 'exon')")
    def slop = task.ext.slop?:""
    verify(!isBlank(slop),"${task.process} : task.ext.slop must be defined (e.g: 1000)")
    def extra_cmd = task.ext.extra_cmd?:""
    verify(!isBlank(extra_cmd),"${task.process} : task.ext.extra_cmd must be defined (e.g: grep -F protein_coding or just cat )") 

	def prefix = task.ext.prefix?:"${gtf.baseName}.${what}s.x${slop}"
"""
hostname 1>&2
mkdir -p TMP


${gtf.name.endsWith(".gz")?"gunzip -c":"cat"} "${gtf}" |\\
    ${extra_cmd} |\\
    awk -F '\t' '/^#/ {next;} (\$3=="${what}") {printf("%s\t%d\t%s\\n",\$1,int(\$4)-1,\$5);}' |\\
    bedtools slop -b ${slop} -g "${fai}" |\\
    LC_ALL=C sort --buffer-size=${task.memory.mega}M -T TMP -t '\t' -k1,1 -k2,2n |\\
    bedtools merge > ${prefix}.bed

cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(tabix version | awk '(NR==1) {print \$NF;}')"
    slop: ${slop}
    type: ${what}
END_VERSIONS
"""

stub:
    def what = task.ext.what?:""
    verify(!isBlank(what),"${task.process} : task.ext.what must be defined (e.g: 'exon')")
    def slop = task.ext.slop?:""
    verify(!isBlank(slop),"${task.process} : task.ext.slop must be defined (e.g: 1000)")
    def extra_cmd = task.ext.extra_cmd?:""
    verify(!isBlank(extra_cmd),"${task.process} : task.ext.extra_cmd must be defined (e.g: grep -F protein_coding or just cat )") 

	def prefix = task.ext.prefix?:"${gtf.baseName}.${what}s.x${slop}"
"""
touch versions.yml "${prefix}.bed"
"""
}