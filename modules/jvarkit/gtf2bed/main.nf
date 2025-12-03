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
process GTF_TO_BED {
	label "process_single"
	tag "${meta.id?:""} "
	afterScript "rm -rf TMP"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	input:
		tuple val(meta1),path(dict)
	        tuple val(meta ),path(gtf)
	output:
		tuple val(meta), path("*.bed"),emit:bed
		tuple val(meta), path("*.bed.gz"),path("*.bed.gz.tbi"),emit:tabix
		path("versions.yml"),emit:versions
	script:
		def args = task.ext.args?:""
        def prefix = task.ext.prefix?:"${meta.id}"
        def jvm =  task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData"
        def awk_expr = task.ext.awk_expr?:"(1==1)"
		def cmd1 = task.ext.cmd1?:"cat"
		def cmd2 = task.ext.cmd2?:"cat"
    """
	hostname 1>&2
    mkdir -p TMP
	
    ${gtf.name.endsWith(".gz")?"gunzip -c ":"cat"} ${gtf} |\\
	${cmd1} |\\
        awk -F '\t' '${awk_expr}' |\\
	jvarkit  ${jvm}  gtf2bed -R "${dict}" ${args} |\\
	${cmd2} |\\
	sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n |\\
	uniq > TMP/jeter.bed
    
    cp TMP/jeter.bed ${prefix}.bed

    bgzip TMP/jeter.bed
    tabix -p bed -f  TMP/jeter.bed.gz

    mv TMP/jeter.bed.gz ${prefix}.bed.gz
    mv TMP/jeter.bed.gz.tbi ${prefix}.bed.gz.tbi


cat << EOF > versions.yml
${task.process}:
	jvarkit: "\$(jvarkit --version)"
EOF
	"""
	
stub: 
def prefix = task.ext.prefix?:"${meta.id}"
"""
mkdir versions.yml ${prefix}.bed.gz ${prefix}.bed.gz.tbi ${prefix}.bed
"""
}
