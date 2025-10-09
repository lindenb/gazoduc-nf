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

process SAMTOOLS_VIEW {
   tag "${meta.id} / ${bam.name}"
   label "process_single"
   afterScript "rm -rf TMP"
   conda "${moduleDir}/../../../conda/bioinfo.01.yml"
   input:
		tuple val(meta1),path(fasta)
		tuple val(meta2),path(fai)
		tuple val(meta4),path(optional_bed)
		tuple val(meta5),path(optional_read_names)
	    tuple val(meta),path(bam),path(bai)
   output:
		tuple val(meta),path("*am"),path("*ai"),optional:true,emit:bam
		path("versions.yml"),emit:versions
   script:
		def prefix = task.ext.prefix?: "${meta.id}.\${MD5}"
		def args1  = task.ext.args1?:""
		def args2  = task.ext.args2?:""
		def args3  = task.ext.args3?:""
		def fmt =  task.ext.format?:"BAM"
		def suffix1 = fmt.toUpperCase().contains("B")?"bam":"cram"
		def suffix2 = fmt.toUpperCase().contains("B")?"bai":"crai"
   """
	hostname 1>&2

	mkdir -p TMP
	samtools view \\
		${args1} ${args2} ${args3} \\
		${fasta?"--reference \"${fasta}\"":""} \\
		--threads ${task.cpus} \\
		${optional_bed?"-L \"${optional_bed}\"":""} \\
		${optional_read_names?"--qname-file \"${optional_read_names}\"":""} \\
		--output-fmt "${fmt}" \\
		-o TMP/jeter.${suffix1} \\
		"${bam}"
		
	SAMPLE=`echo "${bam}" | samtools samples | cut -f 1`
	
	if test "\${SAMPLE}" != "${meta.id}" 
	then
		samtools addreplacerg \\
			${fasta?"--reference \"${fasta}\"":""} \\
			-m overwrite_all \\
			--threads ${task.cpus} \\
			-r "ID:${meta.id}" -r "LB:${meta.id}" -r "SM:${meta.id}" \\
			-o TMP/jeter2.${suffix1} \\
			TMP/jeter.${suffix1}
	
		mv TMP/jeter2.${suffix1} TMP/jeter.${suffix1}
	fi
	
	MD5=`cat TMP/jeter.${suffix1} | md5sum | cut -d ' ' -f1`
	
	samtools index --threads ${task.cpus} TMP/jeter.${suffix1}
	
	mv "TMP/jeter.${suffix1}" "${prefix}.${suffix1}"
	mv "TMP/jeter.${suffix1}.${suffix2}" "${prefix}.${suffix1}.${suffix2}"
	
cat << EOF > versions.yml
"${task.process}":
    samtools:\$(samtools  --version | head -n 1| cut -d ' ' -f2)
EOF
   """
   
   stub:
		def fmt =  task.ext.format?:"BAM"
		def suffix1 = fmt.toUpperCase().contains("B")?"bam":"cram"
		def suffix2 = fmt.toUpperCase().contains("B")?"bai":"crai"
		def prefix = bam.baseName
   """
   touch versions.yml "${prefix}.${suffix1}" "${prefix}.${suffix1}.${suffix2}"
   """
   }
