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

include {isBlank;moduleLoad;getVersionCmd;parseBoolean} from '../../modules/utils/functions.nf'

process GRAPHTYPER_GENOTYPE_01 {
tag "${row.interval?:file(row.bed).name}"
memory "10g"
errorStrategy "retry"
maxRetries 3
cpus 10
afterScript "rm -rf TMP2 TMP"
input:
	val(meta)
	val(graphtyper)
	val(row)
output:
	tuple val(row),path("genotyped.list"),emit:output
	path("version.xml"),emit:version
when:
	true
script:
	if(!row.genomeId) {exit 1,"genomeId missing"}
	if(!row.bams) {exit 1,"bams missing"}
	if(!row.interval && !row.bed) {exit 1,"interval/bed missing"}
	if(row.interval && row.bed) {exit 1,"interval/bed both defined"}

	def genome = params.genomes[row.genomeId]
	def reference = genome.fasta

	def copy_bams = (task.attempt > 1)
	def avg_cov_by_readlen= row.avg_cov_by_readlen?:""
	def arg2 = isBlank(avg_cov_by_readlen)?"":"--avg_cov_by_readlen=${avg_cov_by_readlen}"
"""
hostname 1>&2
${moduleLoad("bcftools bedtools samtools jvarkit")}
mkdir -p TMP
mkdir -p TMP2

if ${row.containsKey("bed")} ; then
	cut -f1,2,3 '${row.bed}' |\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\
		bedtools merge |\
		awk -F '\t' '{printf("%s:%d-%s\\n",\$1,int(\$2)+1,\$3);}' > TMP/input.intervals
	test -s TMP/input.intervals
fi

export TMPDIR=\${PWD}/TMP

if ${copy_bams} ; then
	set -x
	mkdir -p TMP/BAMS
	i=1
	cat "${row.bams}" | while read SAM
	do
		echo "\$i : \${SAM}" 1>&2

		# extract reads in regions
		samtools view --uncompressed -O BAM  ${row.containsKey("bed")?"-M -L ${row.bed}":""} -T "${reference}" "\${SAM}"   ${row.containsKey("bed")?"":"--region=${row.interval}"} |\
			java -jar \${JVARKIT_DIST}/jvarkit.jar samrmdupnames -R "${reference}" --validation-stringency SILENT --samoutputformat BAM --bamcompression 0 |\
			samtools view --reference "${reference}"  -O BAM -o TMP/BAMS/tmp.jeter.bam
		
		samtools index -@ "${task.cpus}" TMP/BAMS/tmp.jeter.bam

		mv -v "TMP/BAMS/tmp.jeter.bam" "TMP/BAMS/tmp.\${i}.bam"
		mv -v "TMP/BAMS/tmp.jeter.bam.bai" "TMP/BAMS/tmp.\${i}.bam.bai"
		echo "TMP/BAMS/tmp.\${i}.bam" >> TMP/bams.list
		i=\$((i+1))
	done

	test -s TMP/bams.list
	
fi

${graphtyper} genotype \
	"${reference}" \
	--output=TMP2 \
	--force_no_copy_reference \
	--force_use_input_ref_for_cram_reading \
	--sams=${copy_bams?"TMP/bams.list":"${row.bams}"} \
	${row.containsKey("bed")?  "--region_file=TMP/input.intervals" : "--region=${row.interval}"} \
	--threads=${task.cpus} \
	${arg2}


rm -rf "TMP2/input_sites"


# bug graphtyper ? https://github.com/DecodeGenetics/graphtyper/issues/136 , reorder samples

find TMP2 -type f -name "*.vcf.gz"  | while read F
do
	bcftools query -l "\${F}" | sort > TMP/ordered.samples.txt
	bcftools view  --threads ${task.cpus} -O b -o "\${F}.bcf"  --samples-file TMP/ordered.samples.txt "\${F}"
	bcftools index --threads ${task.cpus} "\${F}.bcf"

	rm "\${F}" "\${F}.tbi"
done


mv TMP2 OUT
find \${PWD}/OUT -type f -name "*.bcf"  > genotyped.list


##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">run graphtyper</entry>
        <entry key="bed">${row.bed?:""}</entry>
        <entry key="interval">${row.interval?:""}</entry>
        <entry key="bams">${row.bams}</entry>
        <entry key="reference">${reference}</entry>
        <entry key="version">${getVersionCmd("bcftools bedtools")}</entry>
</properties>
EOF
"""
}
