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
include {slurpJsonFile;moduleLoad} from '../../modules/utils/functions.nf'
include {hasFeature;isBlank;backDelete} from './annot.functions.nf'

def TAG="MONDO"
def WHATIZ = "The Mondo Disease Ontology (Mondo) aims to harmonize disease definitions across the world. It is a global community effort to reconcile and curate the very many disease resources within a coherent merged ontology. Original versions of Mondo were constructed entirely automatically and used the IDs of source databases and ontologies. Later, additional manually curated cross-ontology axioms were added, and a native Mondo ID system was used to avoid confusion with source databases."

workflow ANNOTATE_MONDO {
	take:
		genomeId
		bed
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
	def jena_version = "5.0.0-rc1"
	def mondo_version = "2024-03-04"
	def ro = "RO:0004021,RO:0004028,RO:0004001,RO:0004004,RO:0004025,RO:0004020,RO:0004003"
	def whatis="MONDO from ${url}"
"""
hostname 1>&2
mkdir -p TMP/TDB
${moduleLoad("htslib jvarkit bedtools")}
set -x

wget --no-check-certificate -O TMP/jeter.zip "https://archive.apache.org/dist/jena/binaries/apache-jena-${jena_version}.zip"
(cd TMP && unzip jeter.zip 1>&2 && rm -v jeter.zip)


wget --no-check-certificate -O TMP/mondo.owl "https://github.com/monarch-initiative/mondo/releases/download/v${mondo_version}/mondo.owl"


./TMP/apache-jena-${jena_version}/bin/tdbloader --loc=./TMP/TDB TMP/mondo.owl


cat << __EOF__ > TMP/jeter.sparql
prefix owl: <http://www.w3.org/2002/07/owl#>
prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>
prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
prefix mondo: <http://purl.obolibrary.org/obo/mondo#>
prefix RO: <http://purl.obolibrary.org/obo/RO_>
PREFIX obo: <http://purl.obolibrary.org/obo/>
PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>


SELECT DISTINCT
	?gene_name
	?entity_label
	?entity_id
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
  ?entity oboInOwl:id ?entity_id .


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



./TMP/apache-jena-5.0.0-rc1/bin/tdbquery --loc=./TMP/TDB --query=TMP/jeter.sparql  --results=TSV |\
	tail -n +2 |\
	tr -d '"' |\
	tr --complement "|A-Za-z0-9._\t\\n" "_" |\
	tr -s "_" |\
	LC_ALL=C sort -T TMP -t \$'\t' -k1,1 > TMP/output.tsv


java -jar \${JVARKIT_DIST}/jvarkit.jar gtf2bed --columns "gtf.feature,gene_name" -R "${reference}" "${genome.gtf}" |\\
	awk -F '\t' '(\$4=="gene")' |\\
	cut -f1,2,3,5 |\\
	LC_ALL=C sort -T TMP -t '\t' -k4,4 |\\
	uniq > TMP/genes.bed

LC_ALL=C  join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,2.2,2.3' TMP/genes.bed TMP/output.tsv |\\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\\
	bgzip > TMP/${TAG}.bed.gz

tabix --force -p bed TMP/${TAG}.bed.gz

mv TMP/${TAG}.bed.gz ./
mv TMP/${TAG}.bed.gz.tbi ./

echo '##INFO=<ID=${TAG}_DISEASE,Number=.,Type=String,Description="Disease names ${WHATIZ}">' > ${TAG}.header
echo '##INFO=<ID=${TAG}_ID,Number=.,Type=String,Description="Diseases ID ${WHATIZ}">' >> ${TAG}.header
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
"""
cat << __EOF__ > ${TAG}.html
<dl>
<dt>${TAG}</dt>
<dd>${WHATIZ}</dd>
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

bcftools annotate -a "${tabix}" -h "${header}" -c "CHROM,FROM,TO,${TAG}_DISEASE,${TAG}_ID" -O b -o TMP/${TAG}.bcf  --merge-logic '${TAG}_DISEASE:unique,${TAG}_ID:unique' '${row.vcf}'
bcftools index TMP/${TAG}.bcf


cat << EOF > TMP/${TAG}.json
{
"vcf"   : "\${PWD}/OUTPUT/${TAG}.bcf",
"index" : "\${PWD}/OUTPUT/${TAG}.bcf.csi",
"bed"   : "${row.bed}"
}
EOF

###
bcftools query -N -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/${TAG}.count
mv -v TMP/${TAG}.* OUTPUT/
${backDelete(row)}
"""
}
