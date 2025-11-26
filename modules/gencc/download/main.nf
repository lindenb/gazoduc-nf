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

process GENCC_DOWNLOAD {
	tag "${meta.id?:fasta.name}"
	afterScript "rm -rf TMP"
	label "process_single"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	input:
		tuple val(meta),path(gtf)
	output:
		tuple val(meta),path("*.bed.gz"), path("*.bed.gz.tbi"), path("*.header"), emit:bed
		path("versions.yml"),emit:versions
	script:
		def url = task.ext.url?:"https://search.thegencc.org/download/action/submissions-export-tsv"
		def prefix = task.ext.prefix?:"${meta.id}.gencc"
	"""
	hostname 1>&2
	export LC_ALL=C
	mkdir -p TMP
	set -x

	jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP gtf2bed  --columns "gtf.feature,gene_name" "${gtf}" |\\
		awk -F '\t' '\$4=="gene" && \$5!="." && \$5!=""' |\\
		cut -f1,2,3,5 |\\
		LC_ALL=C sort --buffer-size=${task.memory.mega}M -t '\t' -k4,4 -T TMP  |\\
		uniq > TMP/genes.bed


    curl -L "${url}" |\\
		tr -d '"' |\\
        awk -F '\t' '{D=\$4;gsub(/[^A-Za-z0-9:]+/,"_",D);M=\$5;gsub(/[^A-Za-z0-9:]+/,"_",M); if(M=="-") M="."; N=\$5;gsub(/[^A-Za-z0-9:]+/,"_",N); if(N=="-") N=".";printf("%s\t%s\t%s\t%s\\n",\$3,D,M,N);}' |\
		LC_ALL=C sort -T TMP -t '\t' -k1,1 |\\
                uniq > TMP/jeter.b

	test -s TMP/jeter.b


	join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,1.4,2.2,2.3,2.4' TMP/genes.bed TMP/jeter.b |\
		LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n | uniq > TMP/jeter.bed

	bgzip TMP/jeter.bed
        tabix --comment '#' -f -p bed TMP/jeter.bed.gz

	mv TMP/jeter.bed.gz ${prefix}.bed.gz
	mv TMP/jeter.bed.gz.tbi gencc.${prefix}.bed.gz.tbi

    echo '##INFO=<ID=GENCC_MONDO,Number=.,Type=String,Description="Mondo database identifier. The Gene Curation Coalition brings together groups engaged in the evaluation of gene-disease validity with a willingness to share data publicly, to develop consistent terminology for gene curation activities, and to facilitate the consistent assessment of genes that have been reported in association with disease ${url}.">' > ${prefix}.header
    echo '##INFO=<ID=GENCC_DISEASE,Number=.,Type=String,Description="Disease Name, The Gene Curation Coalition brings together groups engaged in the evaluation of gene-disease validity with a willingness to share data publicly, to develop consistent terminology for gene curation activities, and to facilitate the consistent assessment of genes that have been reported in association with disease ${url}.">' >> ${prefix}.header
    echo '##INFO=<ID=GENCC_HPO,Number=.,Type=String,Description="The Human Phenotype Ontology identifier. The Gene Curation Coalition brings together groups engaged in the evaluation of gene-disease validity with a willingness to share data publicly, to develop consistent terminology for gene curation activities, and to facilitate the consistent assessment of genes that have been reported in association with disease ${url}.">' >> ${prefix}.header


cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${url}"
END_VERSIONS
"""

stub:
	def prefix="xxx"
"""
touch versions.yml ${prefix}.bed.gz ${prefix}.bed.gz.tbi ${prefix}.header
"""
}
