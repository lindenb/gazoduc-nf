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
include {slurpJsonFile;moduleLoad} from '../../modules/utils/functions.nf'
include {hasFeature;isBlank;backDelete} from './annot.functions.nf'

def TAG="MONDO"

workflow ANNOTATE_MONDO {
	take:
		genomeId
		vcfs /** json: vcf,index,bed */
	main:

             if(hasFeature("mondo")) {
                        source_ch = DOWNLOAD(genomeId)
						annotate_ch = ANNOTATE(source_ch.bed, source_ch.tbi,source_ch.header,vcfs)
                        out1 = annotate_ch.output
                        out2 = annotate_ch.count
                        out3 = MAKE_DOC(genomeId).output
                        }
                else
                    	{
                        out1 = vcfs
                        out2 = Channel.empty()
                        out3 = Channel.empty()
                        }
	emit:
		output = out1
		count = out2
		doc = out3
}

process DOWNLOAD {
afterScript "rm -rf TMP"
memory "5g"
input:
	val(genomeId)
output:
	path("${TAG}.bed.gz"),emit:bed
	path("${TAG}.bed.gz.tbi"),emit:tbi
	path("${TAG}.header"),emit:header
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
	def url = params.annotations,mondo_owl_url
	def jena_version = "5.0.0"
	def mondo_version = "2024-03-04"
	def ro = "RO:0004021,RO:0004028,RO:0004001,RO:0004004,RO:0004025,RO:0004020,RO:0004003"
	def whatis="MONDO from ${url}"
"""
hostname 1>&2
mkdir -p TMP/TDB
${moduleLoad("htslib jvarkit bedtools")}

wget -O TMP/jeter.zip "https://archive.apache.org/dist/jena/binaries/apache-jena-${jena_version}-rc1.zip"
(cd TMP && unzip jeter.zip && rm jeter.zip)


wget -O TMP/mondo.owl "https://github.com/monarch-initiative/mondo/releases/download/v${mondo_version}/mondo.owl"


./TMP/apache-jena-${jena_version}/bin/tdbloader --loc=./TMP/TDB TMP/mondo.owl


cat << __EOF__ > TMP/jeter.sparql
prefix owl: <http://www.w3.org/2002/07/owl#>
prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>
prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
prefix mondo: <http://purl.obolibrary.org/obo/mondo#>
prefix RO: <http://purl.obolibrary.org/obo/RO_>
PREFIX obo: <http://purl.obolibrary.org/obo/>

SELECT DISTINCT
	?gene_name ?gene
	?property_label ?property
	?entity_label  ?entity
WHERE 
{

  VALUES ?property {
    	RO:0004021
        RO:0004028
        RO:0004001
        RO:0004004
        RO:0004025
        RO:0004020
        RO:0004003 
  }
  
  ?entity rdfs:subClassOf [ rdf:type owl:Restriction ;
  owl:onProperty ?property ;
  owl:someValuesFrom ?gene ] . 
  ?entity rdfs:label ?entity_label .


  ?entity rdfs:subClassOf* obo:MONDO_0000001 .

  
 OPTIONAL {
  ?property rdfs:label ?property_label .
 }

 OPTIONAL {
  ?gene rdfs:label ?gene_name .
 }


FILTER (isIRI(?entity) && STRSTARTS(str(?entity), "http://purl.obolibrary.org/obo/MONDO_"))
FILTER (isIRI(?gene) && regex(str(?gene), "hgnc"))
}

__EOF__


./TMP/apache-jena-${jena_version}/bin/tdbquery --loc=./TMP/TDB --query=TMP/jeter.sparql  --results=TSV > TMP/output.tsv


echo '##INFO=<ID=${TAG},Number=0,Type=Flag,Description="${whatis}">' > ${TAG}.header
"""
}



process MAKE_DOC {
executor "local"
input:
        val(genomeId)
output:
	path("${TAG}.html"),emit:output
script:
	def genome = params.genomes[genomeId]
	def url = genome.rmsk_url
"""
cat << __EOF__ > ${TAG}.html
<dl>
<dt>${TAG}</dt>
<dd>UCSC repeat masker intervals. <a href="${url}">${url}</a></dd>
</dl>
__EOF__
"""
}


process ANNOTATE {
tag "${json.name}"
afterScript "rm -rf TMP"
input:
	path(tabix)
	path(tbi)
	path(header)
	path(json)
output:
	path("OUTPUT/${TAG}.json"),emit:output
	path("OUTPUT/${TAG}.count"),emit:count
script:
	def row = slurpJsonFile(json)
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP OUTPUT

bcftools annotate -a "${tabix}" -h "${header}" -c "CHROM,FROM,TO,${TAG}" -O b -o TMP/${TAG}.bcf '${row.vcf}'
bcftools index TMP/${TAG}.bcf


cat << EOF > TMP/${TAG}.json
{
"vcf"   : "\${PWD}/OUTPUT/${TAG}.bcf",
"index" : "\${PWD}/OUTPUT/${TAG}.bcf.csi",
"bed"   : "${row.bed}"
}
EOF

###
bcftools query -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/${TAG}.count
mv -v TMP/${TAG}.* OUTPUT/
${backDelete(row)}
"""
}
