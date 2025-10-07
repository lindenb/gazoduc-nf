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
process FIND_GVCF_BLOCKS {
	label "process_single"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	afterScript "rm -rf TMP"
	tag "${meta.id?:""} ${optional_bed?optional_bed.name:""}"
	input:
		tuple val(meta),path("GVCFS/*"),path(optional_bed)
	output:
		tuple val(meta ),path("*.bed"),path(optional_bed),emit:bed
		path("versions.yml"),emit:versions
	script:
		def args1 = task.ext.args1?:""
		if(args1.trim().isEmpty()) throw new IllegalArgumentException("${task.process} : args1 must be specified in config file");
		def prefix =  task.ext.prefix?:"${meta.id?:""}.blocks"
	"""
	hostname 1>&2
	mkdir -p TMP
	find GVCFS/ -name "*.vcf.gz" | sort -T TMP -V > TMP/gvcfs.list
	test -s  TMP/gvcfs.list
	touch TMP/jeter.bed

	cat TMP/gvcfs.list | while read F
	do
		bcftools index -s "\${F}" | cut -f1 >> TMP/chroms.txt
	done

	cat TMP/chroms.txt | sort -T TMP -V | uniq | while read C
	do

		jvarkit  -XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP findgvcfsblocks \\
			${args1} \\
			${optional_bed?"--bed ${optional_bed}":""} \\
			--chromosome "\${C}" \\
			-o TMP/jeter.interval_list \\
			TMP/gvcfs.list
		
		awk -F '\t' '/^@/ {next;} {printf("%s\t%d\t%s\\n",\$1,int(\$2)-1,\$3);}' TMP/jeter.interval_list >> TMP/jeter.bed
	
	done


if ${optional_bed?true:false}
then
	bedtools intersect -v \\
		-a <(sort -T TMP -t '\t' -k1,1 -k2,2n ${optional_bed} | bedtools merge) \\
		-b <(sort -T TMP -t '\t' -k1,1 -k2,2n TMP/jeter.bed | bedtools merge) > TMP/not.covered.bed

	if test -s TMP/not.covered.bed
	then
		mv TMP/not.covered.bed ./
	fi
fi

mv TMP/jeter.bed ${prefix}.bed

cat << EOF > versions.yml
${task.process}:
	jvarkit: "\$(jvarkit --version)"
EOF
"""

stub:
def prefix =  task.ext.prefix?:"${meta.id?:""}.blocks"
"""
touch ${prefix}.bed
touch versions.yml
"""
}
