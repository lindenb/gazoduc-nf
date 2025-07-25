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

process DISEASES_DOWNLOAD {
	tag "${meta1.id?:fasta.name}"
	afterScript "rm -rf TMP"
	label "process_single"
	conda "${moduleDir}/../../../../conda/bioinfo.01.yml"
	input:
		tuple val(meta1),path(fasta)
		tuple val(meta2),path(fai)
		tuple val(meta3),path(dict)
		tuple val(meta4),path(gtf),path(gtf_tbi)
	output:
		tuple val(meta1),path("*.bed.gz"), path("*.bed.gz.tbi"), path("*.header"), emit:bed
		path("versions.yml"),emit:versions
		path("doc.md"),emit:doc
	script:
		def url = task.ext.url?:"https://download.jensenlab.org/human_disease_textmining_filtered.tsv"
		def treshold = task.ext.treshold ?:4.5
		def TAG = task.ext.tag?:"DISEASES"
		def WHATIZ = "DISEASES is a weekly updated web resource that integrates evidence on disease-gene associations from automatic text mining, manually curated literature, cancer mutation data, and genome-wide association studies. ${url}. Treshold=${treshold}"

	"""
	hostname 1>&2
	export LC_ALL=C
	mkdir -p TMP
		set -x

	jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP gtf2bed  --columns "gtf.feature,gene_name" -R "${fasta}"  "${gtf}" |\\
		awk -F '\t' '\$4=="gene" && \$5!="." && \$5!=""' |\\
		cut -f1,2,3,5 |\\
		LC_ALL=C sort --buffer-size=${task.memory.mega}M -t '\t' -k4,4 -T TMP  |\\
		uniq > TMP/genes.bed


	 wget  -O - "${url}" |\
		awk -F '\t' '(\$5 > ${treshold}) {D=\$3; gsub(/[\\:]/,"_",D); M=\$4; gsub(/[^A-Za-z0-9]+/,"_",M); printf("%s\t%s\t%s\\n",\$2,D,M);}' |\
		LC_ALL=C sort -T TMP -t '\t' -k1,1 > TMP/jeter.b

	test -s TMP/jeter.b


	join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,1.4,2.2,2.3' TMP/genes.bed TMP/jeter.b |\
		LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n | uniq > TMP/jeter.bed

	bgzip TMP/jeter.bed
        tabix --comment '#' -f -p bed TMP/jeter.bed.gz

	mv TMP/jeter.bed.gz diseases.bed.gz
	mv TMP/jeter.bed.gz.tbi diseases.bed.gz.tbi

	echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="${WHATIZ}">' > diseases.header
	echo '##INFO=<ID=${TAG}_DOID,Number=.,Type=String,Description="${WHATIZ}">' >> diseases.header




cat << 'EOF' > doc.md
# annotations:diseases

> DISEASES (https://diseases.jensenlab.org/Search) is a weekly updated web resource that integrates evidence on disease-gene associations from 
> automatic text mining, manually curated literature, cancer mutation data, and genome-wide association studies.

The treshold limit was "${treshold}".

Tags:

 - `INFO/${TAG}` is the name of the disease
 - `INFO/${TAG}_DOID` is the identifier in the "Disease Ontology" 

> The Disease Ontology has been developed  as a standardized ontology  for human disease with the purpose of providing 
> the biomedical community with consistent, reusable and sustainable descriptions of human disease terms, 
> phenotype characteristics.

EOF

cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${URL}"
END_VERSIONS
	"""
	}

