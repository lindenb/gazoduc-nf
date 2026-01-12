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
process VCF_BED_SET_FILTER {
	label "process_single"
	tag "${meta.id?:""} "
	afterScript "rm -rf TMP"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	input:
            tuple val(meta1),path(bed),path(opt_bed_idx)
            tuple val(meta ),path(vcf)
	output:
		tuple val(meta), path("*.vcf.gz"),emit:vcf
		path("versions.yml"),emit:versions
	script:
		def args1 = task.ext.args1?:""
		def args2 = task.ext.args2?:""
                def args3 = task.ext.args3?:""
                def bed_arg = task.ext.bed_arg?:"--blacklist"
                def prefix = task.ext.prefix?:"${meta.id?:"mantamerger"}."
                def jvm =  task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData"
                def jvarkit = task.ext.jvarkit?:"java ${jvm} -jar \${HOME}/jvarkit.jar"// TODO update when conda released
 """
    hostname 1>&2
    mkdir -p TMP
    bcftools view ${args1} "${vcf}" |\\
    ${jvarkit} vcfbedsetfilter  \\
            ${bed_arg} "${bed}" \\
            ${args2} |\\
    bcftools view ${args3} -o z -o TMP/jeter.vcf.gz
           

     mv TMP/jeter.vcf.gz "${prefix}.vcf.gz"
    

cat << EOF > versions.yml
${task.process}:
	jvarkit: "\$(jvarkit --version)"
EOF
"""
	
stub: 
   def prefix = task.ext.prefix?:"${meta.id?:"mantamerger"}."
"""
touch versions.yml ${prefix}.vcf.gz 
"""
}
