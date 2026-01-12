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
process MANTA_MERGER {
	label "process_single"
	tag "${meta.id?:""} "
	afterScript "rm -rf TMP"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	input:
            tuple val(meta1),path(opt_bed)
            tuple val(meta ),path("VCFS/*")
	output:
		tuple val(meta), path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
		path("versions.yml"),emit:versions
	script:
		def args1 = task.ext.args1?:""
		def args2 = task.ext.args2?:""
                def prefix = task.ext.prefix?:"${meta.id?:"mantamerger"}."
                def jvm =  task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData"
                def jvarkit = task.ext.jvarkit?:"java ${jvm} -jar \${HOME}/jvarkit.jar"// TODO update when conda released
 """
	hostname 1>&2
    mkdir -p TMP
    find VCFS/ -name "*.vcf.gz" > TMP/jeter.list
    test -s TMP/jeter.list

    ${jvarkit} mantamerger  \\
            ${opt_bed?"--bed ${opt_bed}":""} \\
            ${args1} \\
            TMP/jeter.list  > TMP/jeter.vcf

     bcftools view ${args2} --threads ${task.cpus} -O z -o TMP/jeter.vcf.gz
     bcftools index --threads ${task.cpus} --force --tbi TMP/jeter.vcf.gz

     mv TMP/jeter.vcf.gz "${prefix}.vcf.gz"
     mv TMP/jeter.vcf.gz.tbi "${prefix}.vcf.gz.tbi"
    

cat << EOF > versions.yml
${task.process}:
	jvarkit: "\$(jvarkit --version)"
EOF
	"""
	
stub: 
   def prefix = task.ext.prefix?:"${meta.id?:"mantamerger"}."
"""
touch versions.yml ${prefix}.vcf.gz ${prefix}.vcf.gz.tbi
"""
}
