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


process SPLIT_N_VARIANTS {
label "process_single"
tag "${meta.id?:""} ${vcf.name} ${optional_bed?optional_bed.name:""}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(optional_bed)
	tuple val(meta ),path(vcf),path(idx)
output:
	tuple val(meta),path("*.vcf.gz",arity:"0..*"),emit:vcf
	tuple val(meta),path("*.tbi",arity:"0..*"),emit:tbi
	tuple val(meta),path("*.MF"),optional:true,emit:manifest
	path("versions.yml"),emit:versions
script:
	def method = task.ext.method?:""
	if(method.trim().isEmpty()) throw new IllegalArgumentException("method undefined for ${task.process}");
	def has_bed = optional_bed?true:false
	def args1 = task.ext.args1?:""
	def prefix = task.ext.prefix?:"${meta.id}.split${optional_bed?".${meta1.id}":""}"
	def jvm = task.ext.jvm?:"-Xmx${task.memory.giga}g  -XX:-UsePerfData  -Djava.io.tmpdir=TMP"
	def awk_expr = task.ext.awk_expr?:"(1==1)"
	def jvarkit = task.ext.jvarkit?:"java ${jvm} -jar \${HOME}/jvarkit.jar"
"""
hostname 1>&2
mkdir -p TMP/OUT
set -x
if ${has_bed}
then
	bcftools index -s "${vcf}" |\\
		awk -F '\t' '${awk_expr} {printf("%s\t0\t%s\\n",\$1,\$2);}' |\\
		LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter1.bed

		
	${has_bed && optional_bed.endsWith(".gz")?"gunzip -c":"cat"} "${optional_bed}" |\\
		LC_ALL=C sort -S '${task.memory.kilo}'  -T TMP -t '\t' -k1,1 -k2,2n |\\
		bedtools merge > TMP/jeter2.bed

	bedtools intersect \\
		-a  TMP/jeter1.bed \\
		-b TMP/jeter2.bed > TMP/jeter3.bed

	# prevent empty bed
	if test ! -s TMP/jeter3.bed ; then
		echo -e 'K_\t0\t1' > TMP/jeter3.bed
	fi
fi

#cat TMP/jeter3.bed

bcftools view ${args1} \\
	${has_bed?" --regions-file TMP/jeter3.bed":""} \\
	"${vcf}" |\\
	${jvarkit} vcfsplitnvariants \\
	--manifest ${prefix}.MF \\
	${method} \\
	-o TMP/OUT/${prefix}

find TMP/OUT/ -type f 1>&2
find TMP/OUT/ -type f -name "*.vcf.gz" | while read F
do
	bcftools index --threads ${task.cpus} -t --force "\${F}"
done

mv TMP/OUT/*.vcf.gz ./ || true
mv TMP/OUT/*.vcf.gz.tbi ./ || true

cat << EOF > versions.yml
${task.process}:
	jvarkit: "\$(${jvarkit} --version)"
EOF
"""

stub:
"""
mkdir  -p TMP/OUT
for X in 1 2 3 4 5 6 7 8 9 10
do
	touch TMP/OUT/\${X}.vcf.gz
	touch TMP/OUT/\${X}.vcf.gz.tbi
done
mv TMP/OUT/*.vcf.gz ./ || true
mv TMP/OUT/*.vcf.gz.tbi ./ || true
touch versions.yml 
"""
}

