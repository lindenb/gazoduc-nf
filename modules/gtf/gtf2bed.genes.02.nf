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
include {getVersionCmd;moduleLoad;parseBoolean} from '../utils/functions.nf'

/**

produce a bed tabix with chrom/start/name/gene_name/gene_id

*/
process GTF_TO_GENES_BED_02 {
tag "${file(gtf).name}"
afterScript "rm -rf TMP"
memory "1g"
input:
	val(meta)
	val(reference)
	val(gtf)
output:
	path("genes.bed.gz"),emit:bed
	path("genes.bed.gz"),emit:tbi
	path("version.xml"),emit:version
script:
	def with_header = parseBoolean(meta.with_header)
	def slop_arg = meta.slop?:""
	"""
	hostname 1>&2
	${moduleLoad("htslib java jvarkit bedtools")}
	set -o pipefail
	mkdir -p TMP

	gunzip -c "${gtf}" |\
			grep -v "^#" |\
			java -Djava.io.tmpdir=TMP -jar "\${JVARKIT_DIST}/jvarkit.jar" gtf2bed --columns "gtf.feature,gene_type,gene_id,gene_name" |\
			awk -F '\t' '(\$4=="gene")' |\
			cut -f 1-3,5- |\
			java -Djava.io.tmpdir=TMP -jar "\${JVARKIT_DIST}/jvarkit.jar" bedrenamechr --column 1 --convert SKIP -R "${reference}" |\
			${slop_arg.isEmpty()?"":"bedtools slop ${slop_arg} -i - -g \"${reference}.fai\" |"} \
			LC_ALL=C sort -T TMP -k1,1 -k2,2n  |\
			awk -F '\t' 'BEGIN{if(${with_header?1:0}){printf("#chrom\tstart\tend\tgene_type\tgene_id\tgene_name\\n");}} {print}' > TMP/jeter.bed

	bgzip TMP/jeter.bed
	tabix --comment '#' -f -p gff TMP/jeter.bed.gz
	mv TMP/jeter.bed.gz "genes.bed.gz"
	mv TMP/jeter.bed.gz.tbi "genes.bed.gz.tbi"

	
	######################################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Extract BED for Genes in a GTF file</entry>
		<entry key="reference">${reference}</entry>
		<entry key="gtf">${gtf}</entry>
		<entry key="slop_arg">${slop_arg}</entry>
		<entry key="with_header">${with_header}</entry>
		<entry key="versions">${getVersionCmd("jvarkit/jvarkit bedtools bgzip tabix ")}</entry>
	</properties>
	EOF
	"""
	}
