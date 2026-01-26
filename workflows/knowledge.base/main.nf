include {JENA_DOWNLOAD} from '../../modules/jena/download'

workflow {
	meta = [id:"knwoledge"]
	JENA_DOWNLOAD(meta)
	GTF2XML([[id:"gtf"],file("/LAB-DATA/GLiCID/projects/BiRD_resources/species/human/GRCh38/gencode.v38.annotation.gtf.gz")])
	GWAS_CATALOG([id:"meta"])
	DOWNLOAD_ONTOLOGY(
		Channel.of(
		[ [id:"hp"],"https://github.com/obophenotype/human-phenotype-ontology/raw/refs/heads/master/hp-base.owl"],
		[ [id:"doid"],"https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/doid.owl"],
		[ [id:"go"],"https://current.geneontology.org/ontology/go.owl"],
		[ [id:"mondo"] , "https://github.com/monarch-initiative/mondo/releases/download/v2026-01-06/mondo.owl"],
		[ [id:"bto"], "http://purl.obolibrary.org/obo/bto.owl" ],
		[ [id:"uberon"], "http://purl.obolibrary.org/obo/uberon.owl"]
		))
	DOWNLOAD_GOA([id:"goa"])
	TDB_LOAD(JENA_DOWNLOAD.out.tdbloader, 
		GTF2XML.out.rdf
			.mix(GWAS_CATALOG.out.rdf)
			.mix(DOWNLOAD_GOA.out.rdf)
			.mix(DOWNLOAD_ONTOLOGY.out.rdf)
			.map{meta,f->f}.collect().map{f->[[id:"x"],f.sort()]} 
		)
	}


process DOWNLOAD_GOA {

input:
	val(meta)
output:
        tuple val(meta),path("*.gz"),emit:rdf
script:
"""

cat << EOF > goa.rdf
<rdf:RDF
        xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
        xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#"
        xmlns:dc="http://purl.org/dc/elements/1.1/"
        xmlns:u="https://umr1087.univ-nantes.fr/rdf/"
  >
EOF

cat << EOF | awk -F '|' '{printf("<u:evidence rdf:about=\\"https://umr1087.univ-nantes.fr/rdf/evidence/%s\\"><dc:title>%s</dc:title></u:evidence>\\n",\$1,\$2);}'  >> goa.rdf
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
	 awk  '{printf("<u:goa_qualifier rdf=\\"https://umr1087.univ-nantes.fr/rdf/goa/%s\\"><dc:label>%s</dc:label></u:goa_qualifier>\\n",\$1,\$1);}'  >> goa.rdf

cat << 'EOF' > jeter.awk
/^!/ {next;}
(\$4 ~ /^NOT/) {next;}
(\$13!="taxon:9606") {next;}

	{
	G=\$5;
	gsub(/:/,"_",G);
	printf("<u:goa>");
	  printf("<u:gene_name>%s</u:gene_name>",\$3);
	  printf("<u:has_go rdf:resource=\\"http://purl.obolibrary.org/obo/%s\\"/>",G);
	  printf("<u:has_qualifier rdf:resource=\\"https://umr1087.univ-nantes.fr/rdf/goa/%s\\"/>",\$4);
	  printf("<u:has_evidence rdf:resource=\\"https://umr1087.univ-nantes.fr/rdf/evidence/%s\\"/>",\$7);
	printf("</u:goa>\\n");
	}


EOF

wget -O - "https://current.geneontology.org/annotations/goa_human.gaf.gz" | gunzip -c |awk -F '\t' -f jeter.awk >> goa.rdf


echo '</rdf:RDF>' >> goa.rdf
gzip --best goa.rdf
"""
}

process DOWNLOAD_ONTOLOGY {
tag "${u}"
input:
	tuple val(meta),val(u)
output:
	tuple val(meta),path("*.gz"),emit:rdf
script:
"""
wget -O ${meta.id}.rdf "${u}"
gzip "${meta.id}.rdf"
"""
}

process GWAS_CATALOG {
input:
	val(meta)
output:
	 tuple val(meta),path("*.gz"),emit:rdf
script:
"""

cat << 'EOF' > jeter.awk

	{
	if(\$37=="") next;
	gsub(/[^A-Za-z, 0-9\\-;]+/,"_",\$16);
	gsub(/[^A-Za-z, 0-9\\-;]+/,"_",\$17);
	gsub(/[^A-Za-z, 0-9\\-;]+/,"_",\$18);
	gsub(/\\.[0-9]+\$/,"",\$4);
	printf("<u:gwasPeak>");
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
	printf("<u:p-value  rdf:datatype=\\"http://www.w3.org/2001/XMLSchema#double\\">%s</u:p-value>",\$28);
	printf("<u:study>");
		printf("<u:gwasStudy  rdf:about=\\"https://www.ebi.ac.uk/gwas/studies/%s\\">",\$4);
		printf("<dc:title>%s</dc:title>",\$7);
		gsub(/[ ]/,"",\$36);
		N=split(\$36,a,/[,]/);
		for(i=1;i<=N;i++) printf("<u:mappedTrait rdf:resource=\\"%s\\"/>",a[i]);
		printf("</u:gwasStudy>");
	printf("</u:study>");
	printf("</u:gwasPeak>\\n");
	}

EOF

cat << EOF > gwas.rdf
<rdf:RDF
        xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
        xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#"
        xmlns:dc="http://purl.org/dc/elements/1.1/"
        xmlns:u="https://umr1087.univ-nantes.fr/rdf/"
  >
EOF

wget -O jeter.zip "https://www.ebi.ac.uk/gwas/api/search/downloads/associations/v1.0.2?split=false"

unzip -p jeter.zip gwas-catalog-download-associations-alt-full.tsv |\\
	tail -n +2 |\
	head -n 100 |\
	awk -F '\t' -f jeter.awk  >> gwas.rdf

echo "</rdf:RDF>" >> gwas.rdf
gzip gwas.rdf
rm jeter.zip
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
<rdf:RDF
	xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
        xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#"
	xmlns:dc="http://purl.org/dc/elements/1.1/"
        xmlns:u="https://umr1087.univ-nantes.fr/rdf/"
  >
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
	 			printf("<u:build>hg38</u:build>");
	 			printf("<u:chrom>%s</u:chrom>",\$1);
	 			printf("<u:start rdf:datatype=\\"http://www.w3.org/2001/XMLSchema#int\\">%s</u:start>",\$2);
	 			printf("<u:end rdf:datatype=\\"http://www.w3.org/2001/XMLSchema#int\\">%s</u:end>",\$3);
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

		if(\$5!=".")  printf("<u:name>%s</u:name>",\$5);
		if(\$6!=".")  printf("<u:type>%s</u:type>",\$6);
		printf("<u:location>");
			printf("<u:Location>");
	 			printf("<u:build>hg38</u:build>");
	 			printf("<u:chrom>%s</u:chrom>",\$1);
	 			printf("<u:start rdf:datatype=\\"http://www.w3.org/2001/XMLSchema#int\\">%s</u:start>",\$2);
	 			printf("<u:end rdf:datatype=\\"http://www.w3.org/2001/XMLSchema#int\\">%s</u:end>",\$3);
			printf("</u:Location>");
		printf("</u:location>");
		printf("</u:Gene>");
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
