/*

Copyright (c) 2022 Pierre Lindenbaum

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

/**

https://github.com/DecodeGenetics/graphtyper

graphtyper is a graph-based variant caller capable of genotyping population-scale short read data sets. It represents a reference genome and known variants of a genomic region using an acyclic graph structure (a "pangenome reference"), which high-throughput sequence reads are re-aligned to for the purpose of discovering and genotyping SNPs, small indels, and structural variants.

*/
process GRAPHTYPER_DOWNLOAD_01 {
input:
	val(meta)
output:
	path("graphtyper"),emit:executable
	path("version.xml"),emit:version
script:
	def v = (meta.graptyper_version?:"2.7.5")
"""
wget -O graphtyper "https://github.com/DecodeGenetics/graphtyper/releases/download/v${v}/graphtyper"
chmod a+x graphtyper

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">download graphtyper</entry>
        <entry key="version">${v}</entry>
</properties>
EOF
"""
}


process graphTyperGenotype01 {
tag "${file(meta.bed).name} ${file(meta.bams)}"
afterScript "rm -rf results TMP"
label "graphtyper_genotype"
input:
	val(graphtyper)
	val(meta)
output:
	tuple val(meta),path("genotyped.bcf")
script:
	if(!meta.reference) {exit 1,"reference missing"}
	if(!meta.bams) {exit 1,"bams missing"}
	if(!meta.bed) {exit 1,"bed missing"}
"""
hostname 1>&2
module load bcftools/0.0.0
mkdir TMP
export TMPDIR=\${PWD}/TMP

awk -F '\t' '{printf("%s:%d-%d\\n",\$1,int(\$2)+1,\$3);}' "${meta.bed}" > TMP/jeter.regions

${graphtyper} genotype \
	"${meta.reference}" \
	--force_no_copy_reference \
	--force_use_input_ref_for_cram_reading \
	--sams=${meta.bams} \
	--region_file=TMP/jeter.regions \
	--threads=${task.cpus}

find \${PWD}/results/ -type f -name "*.vcf.gz" | grep -v '/input_sites/' > TMP/vcf.list

bcftools concat --file-list TMP/vcf.list \
	--allow-overlaps --remove-duplicates \
	--threads ${task.cpus} -O b -o "genotyped.bcf"

bcftools index --threads ${task.cpus} "genotyped.bcf"
"""
}
