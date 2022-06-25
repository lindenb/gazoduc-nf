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
include {parseBoolean;getKeyValue;getModules;getGencodeGtfUrl;isUrl} from '../utils/functions.nf'

process GTF_TO_GENES_BED_01 {
tag "${file(reference).name}"
afterScript "rm -rf TMP"
memory "1g"
input:
	val(meta)
	val(reference)
	val(jar)
output:
	path("${file(reference).getSimpleName()}.genes.bed${parseBoolean(getKeyValue(meta,"with_tabix","false"))?".gz":""}"),emit:bed
	path("${file(reference).getSimpleName()}.genes.bed.gz.tbi"),emit:tbi,optional:true
	path("version.xml"),emit:version
script:
	def url0 = getKeyValue(meta,"gtf_url",getGencodeGtfUrl(reference))
	def gene_type_regex = getKeyValue(meta,"gene_type_regex",".")
	def gene_id_regex = getKeyValue(meta,"gene_id_regex",".")
	def gene_name_regex = getKeyValue(meta,"gene_name_regex",".")
	def slop_arg = getKeyValue(meta,"slop_arg","")
	def with_header = parseBoolean(getKeyValue(meta,"with_header","false"))
	def with_tabix = parseBoolean(getKeyValue(meta,"with_tabix","false"))
	"""
	hostname 1>&2
	module load ${getModules("htslib java jvarkit bedtools")}
	set -o pipefail
	mkdir TMP
	
	${isUrl(url0)?"wget -O - ":"cat "} "${url0}" |\
			java -Djava.io.tmpdir=TMP -jar "${jar}" --columns "gtf.feature,gene_type,gene_id,gene_name" |\
			awk -F '\t' '(\$4=="gene" && \$5 ~ /${gene_type_regex}/ && \$6 ~ /${gene_id_regex}/ && \$7 ~ /${gene_name_regex}/ )' |\
			cut -f 1-3,5- |\
			java -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/bedrenamechr.jar --column 1 --convert SKIP -R "${reference}" |\
			${slop_arg.isEmpty()?"":"bedtools slop ${slop_arg} -i - -g \"${reference}.fai\" |"} \
			LC_ALL=C sort -T TMP -k1,1 -k2,2n  |\
			awk 'BEGIN{if(${with_header?1:0}){printf("#chrom\tstart\tend\tgene_type\tgene_id\tgene_name\\n");}} {print}' > TMP/jeter.bed

	if [ ! -z "${with_tabix?"Y":""}" ] ; then
		bgzip TMP/jeter.bed
		tabix --comment '#' -p bed TMP/jeter.bed.gz
		mv TMP/jeter.bed.gz "${file(reference).getSimpleName()}.genes.bed.gz"
		mv TMP/jeter.bed.gz.tbi "${file(reference).getSimpleName()}.genes.bed.gz.tbi"
	else
		mv TMP/jeter.bed "${file(reference).getSimpleName()}.genes.bed"
	fi
	
	######################################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Download GTF for Reference</entry>
		<entry key="reference">${reference}</entry>
		<entry key="url0">${url0}</entry>
		<entry key="gene_type_regex">${gene_type_regex}</entry>
		<entry key="gene_id_regex">${gene_id_regex}</entry>
		<entry key="gene_name_regex">${gene_name_regex}</entry>
		<entry key="slop_arg">${slop_arg}</entry>
		<entry key="with_header">${with_header}</entry>
		<entry key="with_tabix">${with_tabix}</entry>
		<entry key="bedrenamechr">\$(java  -jar \${JVARKIT_DIST}/bedrenamechr.jar --version)</entry>
		<entry key="bgzip">\$(bgzip --version | head -n 1)</entry>
		<entry key="tabix">\$(tabix --version | head -n 1)</entry>
	</properties>
	EOF
	"""
	}
