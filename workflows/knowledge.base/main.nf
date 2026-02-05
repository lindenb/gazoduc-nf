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
include { PREPARE_ONE_REFERENCE                   } from '../../subworkflows/samtools/prepare.one.ref'
include { JENA_DOWNLOAD                           } from '../../modules/jena/download'
include { JENA_ARQ as JENA_ARQ_SO                 } from '../../modules/jena/arq'
include { GTF_TO_XML                              } from '../../modules/jvarkit/gtf2xml'
include { XSLTPROC as SPARQLSO2XSL                } from '../../modules/xsltproc'
include { XSLTPROC as GENE2RDF_XSL                } from '../../modules/xsltproc'
include { DOWNLOAD_GTF_OR_GFF3 as DOWNLOAD_GTF    } from '../../modules/gtf/download'
include { LINUX_SPLIT as LINUX_SPLIT_GTF          } from '../../modules/utils/linuxsplit'

def XSD_NS = "http://www.w3.org/2001/XMLSchema"
def U1087_NS = "https://umr1087.univ-nantes.fr/rdf/"
workflow {
	metadata = [id:"knowledge"]
	versions = Channel.empty()
	multiqc =  Channel.empty()
	
	if(params.fasta==null) {
		log.warn("--fasta missing")
		exit -1
		}

	if(params.gtf==null) {
		log.warn("--gtf missing")
		exit -1
		}
	
	PREPARE_ONE_REFERENCE(
			metadata.plus(
				skip_scatter : true
				),
			Channel.fromPath(params.fasta).map{f->[[id:f.baseName],f]}
			)
  	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)

	JENA_DOWNLOAD(metadata)
	versions = versions.mix(JENA_DOWNLOAD.out.versions)
	
	DOWNLOAD_GTF(PREPARE_ONE_REFERENCE.out.dict)
	versions = versions.mix(DOWNLOAD_GTF.out.versions)

	LINUX_SPLIT_GTF(DOWNLOAD_GTF.out.gtf.map{meta,gtf,_tbi->[meta,gtf]})
	versions = versions.mix(LINUX_SPLIT_GTF.out.versions)

	/*
	 * use scatter_bed to find places to cut gtf file
	 */
	ch1 = LINUX_SPLIT_GTF.out.output
		.map{meta,f->f instanceof List?f:[f]}
		.flatMap()
		.map{gtf->[[id:gtf.name.md5()],gtf]}

	
	DOWNLOAD_ONTOLOGY(
		Channel.of(
		[ [id:"hp"],"https://github.com/obophenotype/human-phenotype-ontology/raw/refs/heads/master/hp-base.owl"],
		[ [id:"doid"],"https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/doid.owl"],
		[ [id:"go"],"https://current.geneontology.org/ontology/go.owl"],
		[ [id:"mondo"] , "https://github.com/monarch-initiative/mondo/releases/download/v2026-01-06/mondo.owl"],
		[ [id:"bto"], "http://purl.obolibrary.org/obo/bto.owl" ],
		[ [id:"uberon"], "http://purl.obolibrary.org/obo/uberon.owl"],
		[ [id:"so"], "https://github.com/The-Sequence-Ontology/SO-Ontologies/raw/refs/heads/master/Ontology_Files/so-simple.owl"]
		))
	versions = versions.mix(DOWNLOAD_ONTOLOGY.out.versions)


	JENA_ARQ_SO(
		JENA_DOWNLOAD.out.arq,
		DOWNLOAD_ONTOLOGY.out.rdf.filter{meta,_rdf->meta.id=="so"},
		[[id:"q"],file("${moduleDir}/query.01.sparql")]
		)
	versions = versions.mix(JENA_ARQ_SO.out.versions)
	
	GTF_TO_XML(
		PREPARE_ONE_REFERENCE.out.dict ,
		ch1.take(2).map{meta,gtf->[meta,gtf,[]]}
		)
	versions = versions.mix(GTF_TO_XML.out.versions)


	SPARQLSO2XSL(
		[[id:"rec.biotype"],file("${moduleDir}/sparql.so2xsl.xsl")],
		JENA_ARQ_SO.out.output
		)
	versions = versions.mix(SPARQLSO2XSL.out.versions)

	GENE2RDF_XSL(
		SPARQLSO2XSL.out.xml.map{meta,xsl->[meta,[file("${moduleDir}/gtf2rdf.xsl"),xsl]]}.first(),
		GTF_TO_XML.out.xml
		)
	versions = versions.mix(GENE2RDF_XSL.out.versions)

	DOWNLOAD_GOA([id:"goa"])
	versions = versions.mix(DOWNLOAD_GOA.out.versions)


	DOWNLOAD_NCBI_INFO([id:"ncbi"])
	versions = versions.mix(DOWNLOAD_NCBI_INFO.out.versions)
	/*
	
	XSLTPROC(
		[[id:"xslt"],file("${moduleDir}/../../src/xsl/gtf2rdf.xsl")],
		GTF_TO_XML.out.xml
		)
	versions = versions.mix(XSLTPROC.out.versions)
	*/
}

/*
	GWAS_CATALOG([id:"gwas"])
	DOWNLOAD_ONTOLOGY(
		Channel.of(
		[ [id:"hp"],"https://github.com/obophenotype/human-phenotype-ontology/raw/refs/heads/master/hp-base.owl"],
		[ [id:"doid"],"https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/doid.owl"],
		[ [id:"go"],"https://current.geneontology.org/ontology/go.owl"],
		[ [id:"mondo"] , "https://github.com/monarch-initiative/mondo/releases/download/v2026-01-06/mondo.owl"],
		[ [id:"bto"], "http://purl.obolibrary.org/obo/bto.owl" ],
		[ [id:"uberon"], "http://purl.obolibrary.org/obo/uberon.owl"],
		[ [id:"so"], "https://github.com/The-Sequence-Ontology/SO-Ontologies/raw/refs/heads/master/Ontology_Files/so-simple.owl"]
		))
	
	DOWNLOAD_STRING([id:"stringdb"])
	DOWNLOAD_NCBI_INFO([id:"ncbi"])
	PROTEIN_ATLAS([id:"proteinatlas"])

	TDB_LOAD(JENA_DOWNLOAD.out.tdbloader, 
		GTF2XML.out.rdf
			.mix(DOWNLOAD_NCBI_INFO.out.rdf)
			.mix(GWAS_CATALOG.out.rdf)
			.mix(PROTEIN_ATLAS.out.rdf)
			.mix(DOWNLOAD_GOA.out.rdf)
			.mix(DOWNLOAD_ONTOLOGY.out.rdf)
			.mix(DOWNLOAD_STRING.out.rdf)
			.map{meta,f->f}.collect().map{f->[[id:"x"],f.sort()]} 
		)
	}
*/

process DOWNLOAD_NCBI_INFO {
label "process_single"
input:
        val(meta)
output:
        tuple val(meta),path("*.gz"),emit:rdf
		path("versions.yml"),emit:versions
script:
	def url = task.ext.url?:"https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
"""

cat << 'EOF' > jeter.awk
(NR>1) {
	N=split(\$6,a,/[|]/);
	for(i=1;i<=N;i++) {
		if(a[i] ~ /^Ensembl\\:/) {
			printf("<u:Gene rdf:about=\\"http://rdf.ebi.ac.uk/resource/ensembl/%s\\"/>", substr(a[i],9) );
				printf("  <u:ncbi_gene_id>%s</u:ncbi_gene_id>",\$2);
				printf("  <u:description><![CDATA[%s]]></u:description>",\$12);

			for(j=1;j<=N;j++) {
				if(i==j) continue;
				if(a[j] ~ /^MIM\\:/) {
					printf("<u:omim>%s</u:omim>",substr(a[j],6));
					}
				else if(a[j] ~ /^HGNC\\:/) {
					printf("<u:hgnc>%s</u:hgnc>",substr(a[j],6));
					}
				}
			printf("</u:Gene>\\n");
			break;
			}
		}
	}
EOF

cat << EOF > string.rdf
<?xml version="1.0" encoding="UTF-8" ?>
<rdf:RDF
        xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
        xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#"
        xmlns:dc="http://purl.org/dc/elements/1.1/"
        xmlns:u="${U1087_NS}"
  >
EOF

curl -L "${url}" |\\
	gunzip -c |\\
	awk -F '\t' -f jeter.awk >> ncbi.rdf

echo '</rdf:RDF>' >> ncbi.rdf
gzip ncbi.rdf
rm jeter.awk

echo << EOF > versions.yml
${task.process}:
	ncbi_gene_info: ${url}
EOF
"""

stub:
"""
touch ncbi.rdf versions.yml
gzip ncbi.rdf
"""
}
/*

>>> 9
$1              Gene : ENSG00000000003
$2         Gene name : TSPAN6
$3            Tissue : Breast
$4   IHC tissue name : Breast
$5         Cell type : myoepithelial cells
$6             Level : Not detected
$7       Reliability : Approved
<<< 9

*/
process PROTEIN_ATLAS {
afterScript "rm -rf TMP"
input:
	val(meta)
output:
        tuple val(meta),path("*.gz"),emit:rdf
script:
	def url = "https://www.proteinatlas.org/download/tsv/normal_ihc_data.tsv.zip"
"""
mkdir -p TMP


cat << EOF > TMP/atlas.rdf
<?xml version="1.0" encoding="UTF-8" ?>
<rdf:RDF
        xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
        xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#"
        xmlns:dc="http://purl.org/dc/elements/1.1/"
        xmlns:u="${U1087_NS}"
  >
EOF


cat << EOF | tr ":" "_" | tr "|" "\t" | sort -T TMP -t '\t' -k1,1 > TMP/jeter.a
Adipose tissue|UBERON:0001013
Brown adipose tissue|UBERON:0001013
Adrenal gland|UBERON:0002369
Appendix|UBERON:0013689
Bone marrow|UBERON:0002371
Breast|UBERON:0000310
Breast|UBERON:0001911
Bronchus|UBERON:0002185
Caudate|UBERON:0001873
Cartilage|UBERON:0002418
Cerebellum|UBERON:0002037
Cerebral cortex|UBERON:0000956
Cervix|UBERON:0000002
Choroid plexus|UBERON:0001886
Colon|UBERON:0001155
Duodenum|UBERON:0002114
Endometrium|UBERON:0001295
Epididymis|UBERON:0001301
Esophagus|UBERON:0001043
Fallopian tube|UBERON:0003889
Gallbladder|UBERON:0002110
Heart muscle|UBERON:0000948
Heart muscle|UBERON:0001133
Hippocampus|UBERON:0002305
Hypothalamus|UBERON:0001898
Dorsal raphe|UBERON:0002043
Efferent ducts|UBERON:0006946
Eye|UBERON:0000970
Hair|UBERON:0000329
Kidney|UBERON:0002113
Liver|UBERON:0002107
Lung|UBERON:0002048
Lactating breast|UBERON:0001911
Lymph node|UBERON:0000029
Nasopharynx|UBERON:0001728
Oral mucosa|UBERON:0003343
Ovary|UBERON:0000992
Pancreas|UBERON:0001264
Parathyroid gland|UBERON:0001132
Placenta|UBERON:0001987
Pituitary gland|UBERON:0000007
Prostate|UBERON:0002367
Retina|UBERON:0000966
Rectum|UBERON:0001052
Salivary gland|UBERON:0001044
Seminal vesicle|UBERON:0000998
Skeletal muscle|UBERON:0001134
Skin|UBERON:0000014
Skin|UBERON:0001003
Skin|UBERON:0000014
Small intestine|UBERON:0002108
Smooth muscle|UBERON:0001135
Soft tissue|UBERON:0003104
Sole of foot|UBERON:0008338
Spleen|UBERON:0002106
Stomach|UBERON:0000945
Testis|UBERON:0000473
Thyroid gland|UBERON:0002046
Thymus|UBERON:0002370
Tonsil|UBERON:0002372
Urinary bladder|UBERON:0001255
Substantia nigra|UBERON:0002038
Vagina|UBERON:0000996
N/A|N/A
EOF

curl -o TMP/jeter.zip "${url}" 

unzip -p  TMP/jeter.zip  normal_ihc_data.tsv |\\
	tail -n +2 |\\
	sort -T TMP -t '\t' -k3,3 > TMP/jeter.b


join -t '\t' -1 1 -2 3 -v 2 TMP/jeter.a TMP/jeter.b | head > TMP/jeter.c
cat TMP/jeter.c
test ! -s TMP/jeter.c

join -t '\t' -1 1 -2 3 -o '1.2,2.1,2.6,2.7'  TMP/jeter.a TMP/jeter.b |\\
	awk -F '\t' '{
		if(\$1=="N/A") next;
		printf("<u:ProteinAtlasIHC>");
		printf("<u:has_tissue rdf:resource=\\"http://purl.obolibrary.org/obo/%s\\"/>",\$1);
		printf("<u:has_gene rdf:resource=\\"http://rdf.ebi.ac.uk/resource/ensembl/%s\\"/>",\$2);
		printf("<u:level>%s</u:level>",\$3);
		printf("</u:ProteinAtlasIHC>\\n");
		}' >> TMP/atlas.rdf


echo '</rdf:RDF>' >> TMP/atlas.rdf

gzip TMP/atlas.rdf
mv  TMP/atlas.rdf.gz ./
"""
}


process DOWNLOAD_STRING {
input:
	val(meta)
output:
        tuple val(meta),path("*.gz"),emit:rdf
script:
	def url = "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz"
	def treshold  = task.ext.treshold?:900
"""

cat << EOF > string.rdf
<?xml version="1.0" encoding="UTF-8" ?>
<rdf:RDF
        xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
        xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#"
        xmlns:dc="http://purl.org/dc/elements/1.1/"
        xmlns:u="${U1087_NS}"
  >
EOF

curl "${url}" |\\
	gunzip -c |\\
	awk  '(\$3 > ${treshold} && \$1!=\$2) {
		gsub(/^9606\\./,"",\$1);
		gsub(/^9606\\./,"",\$2); 
		printf("<u:StringInteraction>");
			printf("<u:protein_1 rdf:resource=\\"http://rdf.ebi.ac.uk/resource/ensembl/%s\\"/>",\$1);
			printf("<u:protein_2 rdf:resource=\\"http://rdf.ebi.ac.uk/resource/ensembl/%s\\"/>",\$2);
			printf("<u:treshold  rdf:datatype=\\"${XSD_NS}#double\\">%s</u:treshold>",\$3);
		printf("</u:StringInteraction>\\n");
	}' >> string.rdf

echo '</rdf:RDF>' >> string.rdf

gzip string.rdf
"""
}

process DOWNLOAD_GOA {
label "process_single"
input:
	val(meta)
output:
        tuple val(meta),path("*.gz"),emit:rdf
		path("versions.yml"),emit:versions
script:
	def url=task.ext.url?:"https://current.geneontology.org/annotations/goa_human.gaf.gz"
"""

cat << EOF > goa.rdf
<?xml version="1.0" encoding="UTF-8" ?>
<rdf:RDF
        xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
        xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#"
        xmlns:dc="http://purl.org/dc/elements/1.1/"
		xmlns:owl ="http://www.w3.org/2002/07/owl#"
        xmlns:u="${U1087_NS}"
  >

 <owl:Class rdf:about="${U1087_NS}GOAEvidence">
		<rdfs:label>goa_evidence</rdfs:label>
		<rdfs:comment>Evidence in GOA (like 'Inferred by Curator')</rdfs:comment>
 </owl:Class>
 
 <owl:Class rdf:about="${U1087_NS}GOAQualifier">
		<rdfs:label>goa_qualifier</rdfs:label>
 </owl:Class> 

 <owl:Class rdf:about="${U1087_NS}GOAEntry">
		<rdfs:label>goa_entry</rdfs:label>
		<rdfs:comment>entry in goa (${url})</rdfs:comment>
 </owl:Class> 
 
EOF

cat << EOF | awk -F '|' '{printf("<u:GOAEvidence rdf:about=\\"${U1087_NS}/goa/evidence/%s\\"><dc:title>%s</dc:title><u:weight rdf:datatype=\\"${XSD_NS}#int\\">%s</<u:weight></u:GOAEvidence>\\n",\$1,\$2,(NR<3?.0:\$3));}'  >> goa.rdf
IC|Inferred by Curator 
IBA|Inferred from Biological aspect of Ancestor 
IBD|Inferred from Biological aspect of Descendant 
IDA|Inferred from Direct Assay 
IEA|Inferred from Electronic Annotation 
EXP|Inferred from Experiment 
IEP|Inferred from Expression Pattern 
IGI|Inferred from Genetic Interaction 
IGC|Inferred from Genomic Context 
HDA|Inferred from High Throughput Direct Assay 
HTP|Inferred from High Throughput Experiment 
HEP|Inferred from High Throughput Expression Pattern 
HGI|Inferred from High Throughput Genetic Interaction 
HMP|Inferred from High Throughput Mutant Phenotype 
IKR|Inferred from Key Residues 
IMP|Inferred from Mutant Phenotype 
IPI|Inferred from Physical Interaction 
IRD|Inferred from Rapid Divergence
RCA|Inferred from Reviewed Computational Analysis 
ISA|Inferred from Sequence Alignment 
ISM|Inferred from Sequence Model 
ISS|Inferred from Sequence or structural Similarity 
ISO|Inferred from Sequence Orthology 
ND|No biological Data available 
NAS|Non-traceable Author Statement 
TAS|Traceable Author Statement 
EOF


echo 'acts_upstream_of,acts_upstream_of_negative_effect,acts_upstream_of_or_within,acts_upstream_of_or_within_negative_effect,acts_upstream_of_or_within_positive_effect,acts_upstream_of_positive_effect,colocalizes_with,contributes_to,enables,involved_in,is_active_in,located_in,part_of' |\
	 tr "," "\\n" |\\
	 awk  '{W=0; printf("<u:GOAQualifier rdf=\\"${U1087_NS}/goa/qualifier/%s\\"><dc:title>%s</dc:title><u:weight rdf:datatype=\\"${XSD_NS}#int\\">%s</<u:weight></u:GOAQualifier>\\n",\$1,\$1,W);}'  >> goa.rdf

cat << 'EOF' > jeter.awk
/^!/ {next;}
(\$4 ~ /^NOT/) {next;}
(\$13!="taxon:9606") {next;}

	{
	G=\$5;
	gsub(/:/,"_",G);
	printf("<u:GOAEntry>");
	  printf("<u:gene_name>%s</u:gene_name>",\$3);
	  printf("<u:has_go rdf:resource=\\"http://purl.obolibrary.org/obo/%s\\"/>",G);
	  printf("<u:has_qualifier rdf:resource=\\"${U1087_NS}/goa/qualifier/%s\\"/>",\$4);
	  printf("<u:has_evidence rdf:resource=\\"${U1087_NS}/goa/evidence/%s\\"/>",\$7);

	N=split(\$6,a,/[|]/);
	for(i=1;i<=N;i++) {
		if(a[i] ~ /^PMID\\:/) {
			printf("<u:has_reference rdf:resource=\\"https://pubmed.ncbi.nlm.nih.gov/%s\\"/>",substr(a[i],6));
			}
		}  
	printf("</u:GOAEntry>\\n");
	}

EOF

curl -L "${url}" | gunzip -c |awk -F '\t' -f jeter.awk >> goa.rdf


echo '</rdf:RDF>' >> goa.rdf
gzip --best goa.rdf

cat << EOF > versions.yml
${task.process}:
	goa: ${url}
EOF
"""
stub:
"""
touch goa.rdf versions.yml
gzip goa.rdf
"""
}

process DOWNLOAD_ONTOLOGY {
label "process_single"
tag "${meta.id}"
input:
	tuple val(meta),val(u)
output:
	tuple val(meta),path("*.gz"),emit:rdf
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:"${meta.id}"
"""
mkdir -p TMP
curl -L -o TMP/jeter.rdf "${u}"
gzip TMP/jeter.rdf
mv TMP/jeter.rdf.gz "${prefix}.rdf.gz"
rmdir TMP

cat << EOF > versions.yml
"${task.process}":
	url: ${u}
EOF
"""
stub:
	def prefix = task.ext.prefix?:"${meta.id}"
"""
touch versions.yml "${prefix}.rdf.gz"
"""
}

process GWAS_CATALOG {
input:
	val(meta)
output:
	 tuple val(meta),path("*.gz"),emit:rdf
	 path("versions.yml"),emit:versions
script:
	def url=task.ext.url?"https://www.ebi.ac.uk/gwas/api/search/downloads/associations/v1.0.2?split=false"
	def prefix = task.ext.prefix?:"${meta.id}"
"""

cat << 'EOF' > jeter.awk

	{
	if(\$37=="") next;
	gsub(/[^A-Za-z, 0-9\\-;]+/,"_",\$16);
	gsub(/[^A-Za-z, 0-9\\-;]+/,"_",\$17);
	gsub(/[^A-Za-z, 0-9\\-;]+/,"_",\$18);
	gsub(/\\.[0-9]+\$/,"",\$4);
	printf("<u:GWASCatalogEntry>");
	if(\$16!="" && \$17!="") {
		gsub(/[ ]/,"",\$16);
		N=split(\$16,a,/[,]/);
		for(i=1;i<=N;i++) printf("<u:upstream_gene  rdf:resource=\\"http://rdf.ebi.ac.uk/resource/ensembl/%s\\"/>",a[i]);
		gsub(/[ ]/,"",\$17);
		N=split(\$17,a,/[,]/);
		for(i=1;i<=N;i++) printf("<u:downstream_gene  rdf:resource=\\"http://rdf.ebi.ac.uk/resource/ensembl/%s\\"/>",a[i]);
		}
	else
		{
		gsub(/[ ]/,"",\$18);
		N=split(\$18,a,/[,]/);
		for(i=1;i<=N;i++) {
			printf("<u:upstream_gene  rdf:resource=\\"http://rdf.ebi.ac.uk/resource/ensembl/%s\\"/>",a[i]);
			printf("<u:downstream_gene  rdf:resource=\\"http://rdf.ebi.ac.uk/resource/ensembl/%s\\"/>",a[i]);
			}
		}
	printf("<u:p-value  rdf:datatype=\\"${XSD_NS}#double\\">%s</u:p-value>",\$28);
	printf("<u:study>");
		printf("<u:GWASCatalogStudy  rdf:about=\\"https://www.ebi.ac.uk/gwas/studies/%s\\">",\$4);
		printf("<dc:title>%s</dc:title>",\$7);
		gsub(/[ ]/,"",\$36);
		N=split(\$36,a,/[,]/);
		for(i=1;i<=N;i++) printf("<u:mappedTrait rdf:resource=\\"%s\\"/>",a[i]);
		printf("</u:GWASCatalogStudy>");
	printf("</u:study>");
	printf("</u:GWASCatalogEntry>\\n");
	}

EOF

cat << EOF > gwas.rdf
<?xml version="1.0" encoding="UTF-8" ?>
<rdf:RDF
        xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
        xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#"
        xmlns:dc="http://purl.org/dc/elements/1.1/"
        xmlns:u="${U1087_NS}"
  >

 <owl:Class rdf:about="${U1087_NS}GWASCatalogEntry">
		<rdfs:label>gwas_catalog_entry</rdfs:label>
		<rdfs:comment>Entry in gwas catalog (${url})</rdfs:comment>
 </owl:Class> 

 <owl:Class rdf:about="${U1087_NS}GWASCatalogStudy">
		<rdfs:label>gwas_catalog_study</rdfs:label>
		<rdfs:comment>study in gwas catalog (${url})</rdfs:comment>
 </owl:Class> 

EOF

wget -O jeter.zip "${url}"

unzip -p jeter.zip gwas-catalog-download-associations-alt-full.tsv |\\
	tail -n +2 |\
	head -n 100 |\
	awk -F '\t' -f jeter.awk  >> gwas.rdf

echo "</rdf:RDF>" >> gwas.rdf
gzip gwas.rdf
rm jeter.zip

cat << EOF > versions.yml
"${task.process}":
	url: ${u}
EOF
"""
stub:
	def prefix = task.ext.prefix?:"${meta.id}"
"""
touch versions.yml "${prefix}.rdf.gz"
"""
}

process GTF2XML {
input:
	tuple val(meta),path(gtf)
output:
	tuple val(meta),path("*.gz"),emit:rdf
script:
"""
cat << EOF > gtf.rdf
<?xml version="1.0" encoding="UTF-8" ?>
<rdf:RDF
	xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
        xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#"
	xmlns:dc="http://purl.org/dc/elements/1.1/"
        xmlns:u="${U1087_NS}"
  >

<u:Build rdf:about="${U1087_NS}/build/hg38">
	<u:ucsc_name>hg38</u:ucsc_name>
	<u:ensembl_name>GRCh38</u:ensembl_name>
</u:Build>

EOF

gunzip -c "${gtf}" |\
	awk -F '\t' '(\$3=="gene")' |\
	java -jar \${HOME}/jvarkit.jar gtf2bed -c 'gene_id,gene_name,gene_type'	|\
	awk -F '\t' '{
		if(\$4==".") next;
		gsub(/\\.[0-9]+\$/,"",\$4);
		printf("<u:Gene rdf:about=\\"http://rdf.ebi.ac.uk/resource/ensembl/%s\\">",\$4);
		printf("<u:id>%s</u:id>",\$4);
		if(\$5!=".")  printf("<u:name>%s</u:name>",\$5);
		if(\$6!=".")  printf("<u:type>%s</u:type>",\$6);
		printf("<u:location>");
			printf("<u:Location>");
	 			printf("<u:build rdf:resource=\\"${U1087_NS}/build/hg38\\"/>");
	 			printf("<u:chrom>%s</u:chrom>",\$1);
	 			printf("<u:start rdf:datatype=\\"${XSD_NS}#int\\">%s</u:start>",\$2);
	 			printf("<u:end rdf:datatype=\\"${XSD_NS}#int\\">%s</u:end>",\$3);
			printf("</u:Location>");
		printf("</u:location>");
		printf("</u:Gene>");
	}' >> gtf.rdf


gunzip -c "${gtf}" |\
	awk -F '\t' '(\$3=="transcript")' |\
	java -jar \${HOME}/jvarkit.jar gtf2bed -c 'gene_id,transcript_id,protein_id' |\
	cut -f4,5,6 |\
	sort -T . |\
	uniq |\
	awk -F '\t' '{
		if(\$1==".") next;
		if(\$2=="." && \$3==".") continue;
		gsub(/\\.[0-9]+\$/,"",\$1);
		gsub(/\\.[0-9]+\$/,"",\$2);
		gsub(/\\.[0-9]+\$/,"",\$3);
		if(\$2!=".") {
			printf("<u:Transcript rdf:about=\\"http://rdf.ebi.ac.uk/resource/ensembl/%s\\">",\$2);
			printf("<u:has_gene rdf:resource=\\"http://rdf.ebi.ac.uk/resource/ensembl/%s\\">",\$1);
			printf("<u:Transcript>\\n");

			if(\$3!=".") {
				printf("<u:Protein rdf:about=\\"http://rdf.ebi.ac.uk/resource/ensembl/%s\\">",\$3);
				printf("<u:has_transcript rdf:resource=\\"http://rdf.ebi.ac.uk/resource/ensembl/%s\\">",\$2);
				printf("<u:Protein>\\n");
				}
			}

	}' >> gtf.rdf



cat << EOF >> gtf.rdf
</rdf:RDF>
EOF

gzip gtf.rdf
"""
}

process TDB_LOAD {
input:
	tuple val(meta1),path(tdbloader)
	tuple val(meta2),path(files)

"""
mkdir -p OUT
${tdbloader.toRealPath()} --loc=OUT  ${files}
"""
}
