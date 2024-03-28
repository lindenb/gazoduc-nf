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
nextflow.enable.dsl=2



include {isBlank;moduleLoad;getVersionCmd} from '../../../modules/utils/functions.nf'


workflow GATK3_HC_MINIKIT {
take:
	genomeId
	bams
	each_bed
main:
	version_ch = Channel.empty()
	call_ch = PER_BED(genomeId, bams, each_bed)
	version_ch = version_ch.mix(call_ch.version)
emit:
        output = call_ch.output
        version= version_ch

}


process PER_BED {
tag "${file(bed).name}"
afterScript "rm -rf TMP"
memory {
	switch(task.attempt) {
		case 1: return '20G';
		case 2: return '50G';
		default: return '75G';
		}
	}
cpus 3
errorStrategy 'retry'
maxRetries 3
input:
	val(genomeId)
	path(bams)
	val(bed)
output:
	tuple path("genotyped.bcf"),path("genotyped.bcf.csi"),emit:output
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
"""
hostname 1>&2
module purge
${moduleLoad("java-jdk/8.0.112")}
${moduleLoad("gatk/3.8 bcftools/0.0.0 samtools/0.0.0")}

mkdir -p TMP/BAMS TMP/TMP

java -version 1>&2

cp -v "${moduleDir}/Minikit3.java" TMP/TMP/Minikit.java
javac -d TMP/TMP -cp ${params.gatkjar} -sourcepath 'TMP/TMP' 'TMP/TMP/Minikit.java'

cat << EOF > TMP/TMP/log4j.properties
# Set root logger level to DEBUG and its only appender to A1.
log4j.rootLogger=ALL, A1

# A1 is set to be a ConsoleAppender.
log4j.appender.A1=org.apache.log4j.ConsoleAppender

# A1 uses PatternLayout.
log4j.appender.A1.layout=org.apache.log4j.PatternLayout
log4j.appender.A1.layout.ConversionPattern=%-4r [%t] %-5p %c %x - %m%n
EOF
jar cvf TMP/minikit.jar -C TMP/TMP .


awk -F '\t' 'BEGIN {printf("FLAGS=-M -L ${bed} -F 3844 ${params.mapq && ((params.mapq as int)>0)?"-q ${params.mapq}":""} --write-index -T ${reference} -O BAM\\n");} {printf("TMP/BAMS/file%d.bam: %s\\n\tsamtools view \$(FLAGS) -o \\"\$@##idx##\$@.bai\\" \$< && touch -c \$(addsuffix .bai,\$@)\\n",NR,\$0);} END {printf("TMP/bams.list: \$(addprefix TMP/BAMS/file,\$(addsuffix .bam,");for(i=1;i<=NR;i++) printf(" %d",i);printf("))\\n");for(i=1;i<=NR;i++) printf("\techo TMP/BAMS/file%d.bam >> \$@\\n",i);}' "${bams}" > TMP/jeter.mk

make -j ${task.cpus} -f TMP/jeter.mk TMP/bams.list

test -s TMP/bams.list


java -XX:ParallelGCThreads=1 -Xmx${task.memory.giga}g -Dsamjdk.compression_level=1 -Djava.io.tmpdir=TMP -cp ${params.gatkjar}:TMP/minikit.jar Minikit \
                -I "TMP/bams.list" \\
                -L "${bed}" \\
		--attempt "${task.attempt}" \\
		--threads "${task.cpus}" \\
                ${params.mapq && ((params.mapq as int)>0)?"--mapq ${params.mapq}":""} \\
                -o TMP/selection.vcf.gz \\
                -R "${reference}" \\
                ${isBlank(genome.dbsnp)?"":"--dbsnp ${genome.dbsnp}"} 1>&2


bcftools view --threads ${task.cpus} --compression-level 9 -O b -o genotyped.bcf TMP/selection.vcf.gz 
bcftools index -f --threads ${task.cpus} genotyped.bcf

#####
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">call bed</entry>
        <entry key="bed">${bed}</entry>
        <entry key="bams">${bams}</entry>
        <entry key="mapq">${params.mapq?:""}</entry>
        <entry key="dbsnp">${genome.dbsnp}</entry>
</properties>
EOF
"""
}
