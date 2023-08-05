/*

Copyright (c) 2023 Pierre Lindenbaum

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

include {getVersionCmd;isBlank;moduleLoad} from '../../modules/utils/functions.nf'
include {CONCAT_FILES_01} from '../../modules/utils/concat.files.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {DOWNLOAD_GNOMAD_SV_01} from '../gnomad/download_gnomad_sv.01.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'

boolean hasFeature(String key) {
	String k = "with_"+key;
	if(!params.annotations.containsKey(k)) {
		log.warn("undefined params.annotations.${k}")
		return false;
		}
	return (params.annotations[k] as boolean)
	}


workflow BCFTOOLS_ANNOTATE_SOURCES {
	take:
		meta
		genomeId
		bed /** limit to that bed or NO_FILE */
	main:
		version_ch = Channel.empty()
		annot_ch  = Channel.empty()

		if(bed.name.equals("NO_FILE")) {
			sorted_bed = bed
			}
		else
			{
			sortbed_ch = SORT_BED([:], genomeId, bed)
			version_ch = version_ch.mix(sortbed_ch.version)
			sorted_bed = sortbed_ch.bed
			}
		


		if (hasFeature("greendb") && params.genomes[genomeId].containsKey("greendb_url")) {
			db1_ch = DOWNLOAD_GREENDB([:], genomeId, sorted_bed)
			version_ch = version_ch.mix(db1_ch.version)
			annot_ch = annot_ch.mix(db1_ch.output)
			}

		if(hasFeature("go") && params.annotations.containsKey("goTerms") && params.genomes[genomeId].containsKey("gff3")) {
			db6_ch = DOWNLOAD_GO([:], genomeId)
			version_ch = version_ch.mix(db6_ch.version)
			annot_ch = annot_ch.mix(db6_ch.output)
			}

		if(hasFeature("remap") && params.genomes[genomeId].containsKey("remap_url")) {
			db7_ch = DOWNLOAD_REMAP([:], genomeId)
			version_ch = version_ch.mix(db7_ch.version)
			annot_ch = annot_ch.mix(db7_ch.output)
			}


		if(hasFeature("simpe_repeats") &&  params.genomes[genomeId].containsKey("simple_repeats_url")) {
			db8_ch = DOWNLOAD_SIMPLE_REPEATS([:], genomeId)
			version_ch = version_ch.mix(db8_ch.version)
			annot_ch = annot_ch.mix(db8_ch.output)
			}

		if(hasFeature("rmsk") && params.genomes[genomeId].containsKey("rmsk_url")) {
			db9_ch = DOWNLOAD_RMSK([:], genomeId)
			version_ch = version_ch.mix(db9_ch.version)
			annot_ch = annot_ch.mix(db9_ch.output)
			}
		
		if(hasFeature("gnomad_sv") && params.genomes[genomeId].containsKey("ucsc_name") &&  (params.genomes[genomeId].ucsc_name.equals("hg19") || params.genomes[genomeId].ucsc_name.equals("hg38"))) {
			svch = DOWNLOAD_GNOMAD_SV_01([:], genomeId)
			version_ch = version_ch.mix(svch.version)

			db10_ch = FREQUENT_GNOMADSV([:], svch.bed)
			version_ch = version_ch.mix(db10_ch.version)
			annot_ch = annot_ch.mix(db10_ch.output)
			}

		if(hasFeature("vista") && params.genomes[genomeId].containsKey("vista_enhancers_url") ) {
			db11_ch = DOWNLOAD_VISTA_ENHANCERS([:], genomeId)
			version_ch = version_ch.mix(db11_ch.version)
			annot_ch = annot_ch.mix(db11_ch.output)
			}

		
		if(hasFeature("regfeatures") && params.genomes[genomeId].containsKey("ensembl_regulatory_gff_url")) {
			db12_ch = DOWNLOAD_REGFEATURES([:], genomeId)
			version_ch = version_ch.mix(db12_ch.version)
			annot_ch = annot_ch.mix(db12_ch.output)
			}


		if(hasFeature("encode_ccre") && params.genomes[genomeId].containsKey("encode_ccre_url")) {
			db13_ch = DOWNLOAD_ENCODE_CCRE([:], genomeId, bed)
			version_ch = version_ch.mix(db13_ch.version)
			annot_ch = annot_ch.mix(db13_ch.output)
			}

		if(hasFeature("clinvar") && params.genomes[genomeId].containsKey("clinvar_vcf_url")) {
			db13_ch = DOWNLOAD_CLINVAR([:], genomeId)
			version_ch = version_ch.mix(db13_ch.version)
			annot_ch = annot_ch.mix(db13_ch.vcf)
			}



		if (params.genomes[genomeId].containsKey("gtf")) {
	
			bed4genes_ch = BED_FOR_GENES([:], genomeId , sorted_bed)
	                version_ch = version_ch.mix(bed4genes_ch.version)

			if(hasFeature("bhfucl")) {
				db2_ch = DOWNLOAD_BHFUCL([:], bed4genes_ch.bed)
				version_ch = version_ch.mix(db2_ch.version)
				annot_ch = annot_ch.mix(db2_ch.output)
				}

			if(hasFeature("diseases")) {
				db3_ch = DOWNLOAD_DISEASES([:], bed4genes_ch.bed)
				version_ch = version_ch.mix(db3_ch.version)
				annot_ch = annot_ch.mix(db3_ch.output)
				}

			if(hasFeature("gencc")) {
				db4_ch = DOWNLOAD_GENCC([:], bed4genes_ch.bed)
				version_ch = version_ch.mix(db4_ch.version)
				annot_ch = annot_ch.mix(db4_ch.output)
				}


			if(hasFeature("organizer")) {
				db5_ch = DOWNLOAD_ORGANIZER([:], bed4genes_ch.bed)
				version_ch = version_ch.mix(db5_ch.version)
				annot_ch = annot_ch.mix(db5_ch.output)
				}
			}

		annot_ch = annot_ch.map{T->[
			name: T[0],
                        tabix : T[1],
                        header: T[2],
                        columns: T[3],
                        merge_logic: T[4],
                        minoverlap : "."
                        ]}


		version_ch = MERGE_VERSION("annotations.sources", version_ch.collect())		
	emit:
		version = version_ch
		annotations_ch = annot_ch
	}


process SORT_BED {
tag "${bed}"
input:
	val(meta)
	val(genomeId)
	path(bed)
output:
	path("sorted.bed"),emit:bed
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
"""
hostname 1>&2
${moduleLoad("bedtools jvarkit")}
set -o pipefail

${bed.name.endsWith(".gz")?"gunzip -c ":"cat "} "${bed}" |\
	grep -vE "^(#|browser|track)" |\
	cut -f1,2,3 |\
	java -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${reference}" --column 1 --convert SKIP  |\
	sort -T . -t '\t' -k1,1 -k2,2n |\
	bedtools merge > sorted.bed

######
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
	<entry key="description">sort bed</entry>
        <entry key="reference">${reference}</entry>
        <entry key="bed">${bed}</entry>
        <entry key="versions">${getVersionCmd("bedtools jvarkit/bedrenamechr")}</entry>
</properties>
"""
}


/* gene symbols as bed file sorted on gene symbol */
process BED_FOR_GENES {
afterScript "rm -rf TMP"
input:
	val(meta)
	val(genomeId)
	path(sortedbed)
output:
       	path("genes.symbols.${genomeId}.bed"),emit:bed
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
	def gtf = genome.gtf
"""
hostname 1>&2
${moduleLoad("bedtools jvarkit")}
set -o pipefail
mkdir -p TMP


java -jar \${JVARKIT_DIST}/gtf2bed.jar --columns "gtf.feature,gene_name" -R "${reference}" "${gtf}" |\
awk -F '\t' '(\$4=="gene")' |\
cut -f1,2,3,5 |\
${sortedbed.name.equals("NO_FILE")?"":" sort -T TMP -t '\t' -k1,1 -k2,2n | bedtools intersect -u -a - -b '${sortedbed}' | "} \
sort -T TMP -t '\t' -k4,4 | uniq > TMP/genes.symbols.${genomeId}.bed

mv "TMP/genes.symbols.${genomeId}.bed" ./

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
	<entry key="description">convert gtf to bed</entry>
        <entry key="gtf">${gtf}</entry>
        <entry key="bed">${sortedbed}</entry>
        <entry key="versions">${getVersionCmd("jvarkit/gtf2bed awk bedtools")}</entry>
</properties>
"""
}

/* Cardiovascular Gene Ontology Annotation Initiative **/
process DOWNLOAD_BHFUCL {
afterScript "rm -rf TMP"
input:
	val(meta)
      	path(bed)
output:
	tuple val("BHFUCL"),path("bhfucl.bed.gz"),path("bhfucl.header"),val("CHROM,FROM,TO,BHFUCL"),val(""),emit:output
        path("version.xml"),emit:version
when:
script:
	def url =  "http://ftp.ebi.ac.uk/pub/databases/GO/goa/bhf-ucl/gene_association.goa_bhf-ucl.gz"
	def whatis= "Cardiovascular Gene Ontology Annotation Initiative https://www.ebi.ac.uk/GOA/CVI"
"""
hostname 1>&2
${moduleLoad("htslib")}
set -o pipefail
mkdir -p TMP

wget -O - "${url}" |\
	gunzip -c |\
	awk -F '\t' '/#/{next} (\$4 ~ /^NOT/) {next;} {print \$3;}' |\
	uniq | sort -T TMP | uniq |
	sed 's/\$/\t1/' > TMP/genes.txt

join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,2.2' "${bed}" TMP/genes.txt |\
sort -T TMP -t '\t' -k1,1 -k2,2n | bgzip > TMP/bhfucl.bed.gz

tabix -p bed TMP/bhfucl.bed.gz

mv TMP/bhfucl.bed.gz ./
mv TMP/bhfucl.bed.gz.tbi ./

echo '##INFO=<ID=BHFUCL,Number=0,Type=Flag,Description="${whatis}">' > bhfucl.header

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">${whatis}</entry>
        <entry key="versions">${getVersionCmd("tabix wget")}</entry>
        <entry key="url"><a>${url}</a></entry>
</properties>
EOF
"""
}


/* DISEASES is a weekly updated web resource that integrates evidence on disease-gene associations from automatic text mining, manually curated literature, cancer mutation data, and genome-wide association studies. https://diseases.jensenlab.org/ */
process DOWNLOAD_DISEASES {
afterScript "rm -rf TMP"
input:
	val(meta)
      	path(bed)
output:
	tuple val("DISEASES"),path("diseases.bed.gz"),path("diseases.header"),val("CHROM,FROM,TO,DISEASES"),val("DISEASES:unique"),emit:output
	path("version.xml"),emit:version
script:
	def url = "https://download.jensenlab.org/human_disease_knowledge_filtered.tsv"
	def whatis="DISEASES is a weekly updated web resource that integrates evidence on disease-gene associations from automatic text mining, manually curated literature, cancer mutation data, and genome-wide association studies. https://diseases.jensenlab.org/"

if(!isBlank(url))
"""
hostname 1>&2
${moduleLoad("htslib")}
set -o pipefail
mkdir -p TMP

wget -O - "https://download.jensenlab.org/human_disease_knowledge_filtered.tsv" |\
	cut -f 2,4 | tr -c "[a-zA-Z_0-9\\-\\n\\t.]" "_" | tr -s "_"  |\
	sort -T TMP  -t '\t' -k1,1  | uniq > TMP/org.tsv

join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,2.2' "${bed}" TMP/org.tsv |\
sort -T TMP -t '\t' -k1,1 -k2,2n | bgzip > TMP/diseases.bed.gz

tabix -p bed TMP/diseases.bed.gz

mv TMP/diseases.bed.gz ./
mv TMP/diseases.bed.gz.tbi ./

echo '##INFO=<ID=DISEASES,Number=.,Type=String,Description="${whatis}">' > diseases.header


#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">${whatis}</entry>
        <entry key="url"><a>${url}</a></entry>
        <entry key="versions">${getVersionCmd("tabix wget")}</entry>
</properties>
EOF

"""
else
"""
touch diseases.tsv
echo "<properties/>" > version.xml
"""
}


/** gencc The Gene Curation Coalition brings together groups engaged in the evaluation of gene-disease validity with a willingness to share data publicly https://thegencc.org/about.html*/
process DOWNLOAD_GENCC {
afterScript "rm -rf TMP"
input:
	val(meta)
      	path(bed)
output:
	tuple val("GENCC"),path("gencc.bed.gz"),path("gencc.header"),val("CHROM,FROM,TO,GENCC"),val("GENCC:unique"),emit:output
        path("version.xml"),emit:version
script:
	def url= "https://search.thegencc.org/download/action/submissions-export-tsv"
	def whatis = "gencc The Gene Curation Coalition brings together groups engaged in the evaluation of gene-disease validity with a willingness to share data publicly https://thegencc.org/about.html ${url} "

if(!isBlank(url))
"""
hostname 1>&2
${moduleLoad("htslib bedtools")}
set -o pipefail
mkdir -p TMP

wget --no-check-certificate -q  -O - "${url}" |\
	awk -F '\t' '(NR==1) {C=-1;D=-1;for(i=1;i<=NF;i++) {if(\$i=="\\"gene_symbol\\"") C=i;if(\$i=="\\"disease_title\\"") D=i;}next;} {if(C<1 ||  D<1) next; G=\$C;H=\$D; gsub(/[^A-Za-z0-9\\.\\-]+/,"_",G);gsub(/[^A-Za-z0-9\\.\\-]+/,"_",H);  printf("%s\t%s\\n",G,H);}'  |\
	sort -T TMP -t '\t' -k1,1  | uniq > TMP/org.tsv

join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,2.2' "${bed}" TMP/org.tsv |\
sort -T TMP -t '\t' -k1,1 -k2,2n | bgzip > TMP/gencc.bed.gz

tabix -p bed -f TMP/gencc.bed.gz

mv TMP/gencc.bed.gz ./
mv TMP/gencc.bed.gz.tbi ./

echo '##INFO=<ID=GENCC,Number=.,Type=String,Description="${whatis}">' > gencc.header


#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">${whatis}</entry>
        <entry key="url"><a>${url}</a></entry>
        <entry key="versions">${getVersionCmd("tabix wget awk")}</entry>
</properties>
EOF

"""
else
"""
touch gencc.tsv
echo "<properties/>"
"""
}



/** GREEN-DB is a comprehensive collection of 2.4 million regulatory elements in the human genome collected from previously published databases, high-throughput screenings and functional studies.  */
process DOWNLOAD_GREENDB {
afterScript "rm -rf TMP"
input:
	val(meta)
	val(genomeId)
	path(sortedbed)
output:
	tuple val("GREENDB"),path("greendb.bed.gz"),path("greendb.header"),val("CHROM,FROM,TO,GREENDB"),val("GREENDB:unique"),emit:output
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
	def url = genome.greendb_url
	def whatis = "Data from greendb https://zenodo.org/record/5636209.  A comprehensive collection of 2.4 million regulatory elements in the human genome collected from previously published databases, high-throughput screenings and functional studies."
"""
hostname 1>&2
${moduleLoad("htslib bedtools jvarkit")}
set -o pipefail
mkdir -p TMP

wget -O - "${url}" |\
	gunzip -c |\
	grep -v '^chromosome' |\
	cut -f1-3,5 |\
	java -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${reference}" --column 1 --convert SKIP  |\
	${sortedbed.name.equals("NO_FILE")?"":" sort -T . -t '\t' -k1,1 -k2,2n | bedtools intersect -u -a - -b '${sortedbed}' | "} \
	LC_ALL=C sort -T . -t '\t' -k1,1 -k2,2n |\
	uniq |\
	bgzip > TMP/greendb.bed.gz


tabix -p bed TMP/greendb.bed.gz

mv TMP/greendb.bed.gz ./
mv TMP/greendb.bed.gz.tbi ./

echo '##INFO=<ID=GREENDB,Number=.,Type=String,Description="${whatis}. ${url}.">' > greendb.header

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
	<entry key="description">${whatis}</entry>
        <entry key="url"><a>${url}</a></entry>
        <entry key="versions">${getVersionCmd("tabix wget awk")}</entry>
</properties>
EOF
"""
}


process DOWNLOAD_ORGANIZER {
afterScript "rm -rf TMP"
input:
	val(meta)
      	path(bed)
output:
	tuple val("ORGANIZER"),path("organizer.bed.gz"),path("organizer.header"),val("CHROM,FROM,TO,ORGANIZER"),val("ORGANIZER:unique"),emit:output
	path("version.xml"),emit:version
script:

	def url = "http://geneorganizer.huji.ac.il/assets/Uploads/diseases-organs.v12-1.csv.gz"
	def whatis = "Location from Gene Organizer ${url}"

if(!isBlank(url))
"""
hostname 1>&2
${moduleLoad("htslib")}
set -o pipefail
mkdir -p TMP

wget -O - "${url}" | gunzip -c |\
cut -d, -f1,3 | uniq |\
tr -d '"' | tr "," "\t" |\
tr " " "_" | sort -t '\t' -T . -k1,1 > TMP/org.tsv

join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,2.2' "${bed}" TMP/org.tsv |\
sort -T TMP -t '\t' -k1,1 -k2,2n | bgzip > TMP/organizer.bed.gz

tabix -p bed TMP/organizer.bed.gz

mv TMP/organizer.bed.gz ./
mv TMP/organizer.bed.gz.tbi ./

echo '##INFO=<ID=ORGANIZER,Number=.,Type=String,Description="${whatis}">' > organizer.header

echo "\${PWD}/organizer.bed.gz\tCHROM,FROM,TO,ORGANIZER\t\${PWD}/organizer.header" > organizer.tsv

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">${whatis}</entry>
        <entry key="url"><a>${url}</a></entry>
</properties>
EOF
"""
}

process DOWNLOAD_GO {
afterScript "rm -rf TMP"
memory "10g"
input:
	val(meta)
	val(genomeId)
output:
	tuple val("GO"),path("go.bed.gz"),path("go.header"),val("CHROM,FROM,TO,GO"),val(""),emit:output
	path("go.terms.tsv"),emit:goterms_tsv
        path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def gff3 = genome.gff3
	def url0 = "http://purl.obolibrary.org/obo/go/go-basic.obo"
	def url1 = "http://geneontology.org/gene-associations/goa_human.gaf.gz"
	def whatis = "GeneOntology for terms ${meta.goTerms}"

"""
hostname 1>&2
set -o pipefail
${moduleLoad("jvarkit bedtools htslib")}

mkdir -p TMP

wget -O "TMP/go-basic.obo" "${url0}"
wget -O "TMP/goa_human.gaf.gz" "${url1}"

echo "${params.annotations.goTerms}" | tr " ,|;" "\\n" | sort |uniq | grep -v '^\$' > TMP/jeter.terms
test -s TMP/jeter.terms


java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar ${JVARKIT_DIST}/jvarkit.jar goutils \
	--action dump_table \
	--accession-file "TMP/jeter.terms" \
	-go TMP/go-basic.obo  | grep -F -f TMP/jeter.terms > TMP/go.terms.tsv

gunzip -c "${gff3}" | awk -F '\t' '(\$1 ~ /^#/ || \$3=="gene")' | \
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar ${JVARKIT_DIST}/jvarkit.jar goutils \
	--action gff3 \
	--accession-file "TMP/jeter.terms" \
	-go TMP/go-basic.obo \
	-goa TMP/goa_human.gaf.gz  > TMP/jeter.gff

awk '/^#/ {next} {printf("%s\t%d\t%s\\n",\$1,int(\$4)-1,\$5);}' TMP/jeter.gff |\
	sort -T TMP -t '\t' -k1,1 -k2,2n |\
	bedtools merge |\
	sed 's/\$/\t1/' |\
	bgzip > TMP/go.bed.gz

tabix -p bed TMP/go.bed.gz

mv TMP/go.bed.gz ./
mv TMP/go.bed.gz.tbi ./
mv TMP/go.terms.tsv ./


echo -n '##INFO=<ID=GO,Number=0,Type=Flag,Description="${whatis}. Terms: ' > go.header

cat go.terms.tsv |tr  '"' "'" | tr "\t" " " | paste -sd ' ' | tr -d '\\n' >> go.header

echo '">' >> go.header

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">${whatis}</entry>
        <entry key="obo.url"><a>${url0}</a></entry>
        <entry key="annot.url"><a>${url1}</a></entry>
</properties>
EOF

"""
}


/* ReMap is a large scale integrative analysis of DNA-binding experiments for Homo sapiens, Mus musculus, Drosophila melanogaster and Arabidopsis thaliana transcriptional regulators. The catalogues are the results of the manual curation of ChIP-seq, ChIP-exo, DAP-seq from public sources (GEO, ENCODE, ENA).  */
process DOWNLOAD_REMAP {
afterScript "rm -rf TMP"
input:
	val(meta)
	val(genomeId)
output:
	tuple val("REMAP"),path("remap.bed.gz"),path("remap.header"),val("CHROM,FROM,TO,REMAP"),val(""),emit:output
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def url = genome.remap_url
	def whatis="ReMap is a large scale integrative analysis of DNA-binding experiments for Homo sapiens, Mus musculus, Drosophila melanogaster and Arabidopsis thaliana transcriptional regulators. The catalogues are the results of the manual curation of ChIP-seq, ChIP-exo, DAP-seq from public sources (GEO, ENCODE, ENA). ${url}"
	def reference = genome.fasta

"""
hostname 1>&2
${moduleLoad("htslib jvarkit bedtools")}
mkdir -p TMP

set -o pipefail
wget -O - "${url}" |\
	gunzip -c |\
	cut -f1,2,3 |\
	sed 's/\$/\t1/' |\
	java   -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${reference}" --column 1 --convert SKIP  |\
		sort -t '\t' -k1,1 -k2,2n -T TMP |\
		uniq |\
		bgzip > TMP/remap.bed.gz && \
	tabix -p bed -f TMP/remap.bed.gz


mv TMP/remap.bed.gz ./
mv TMP/remap.bed.gz.tbi ./


echo '##INFO=<ID=REMAP,Number=0,Type=Flag,Description="${whatis}">' > remap.header

echo "\${PWD}/remap.bed.gz\tCHROM,FROM,TO,REMAP\t\${PWD}/remap.header" > remap.tsv

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">${whatis}</entry>
        <entry key="url"><a>${url}</a></entry>
</properties>
EOF
"""
}


process DOWNLOAD_SIMPLE_REPEATS {
afterScript "rm -rf TMP"
input:
	val(meta)
	val(genomeId)
output:
	tuple val("SREPEAT"),path("simpleRepeats.bed.gz"),path("simpleRepeats.header"),val("CHROM,FROM,TO,SREPEAT"),val(""),emit:output
	path("version.xml"),emit:version
script:
        def genome = params.genomes[genomeId]
	def reference = genome.fasta
	def url = genome.simple_repeats_url
	def whatis = "Simple repeats from ${url}"

"""
hostname 1>&2
${moduleLoad("htslib jvarkit bedtools")}
mkdir -p TMP

set -o pipefail
wget -O - "${url}" |\
	gunzip -c |\
	cut -f2-4 |\
	sed 's/\$/\t1/' |\
	java -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${reference}" --column 1 --convert SKIP  |\
		sort -t '\t' -k1,1 -k2,2n -T TMP |\
		uniq |\
		bgzip > TMP/simpleRepeats.bed.gz 
	
tabix -p bed -f TMP/simpleRepeats.bed.gz


echo '##INFO=<ID=SREPEAT,Number=0,Type=Flag,Description="${whatis}">' > simpleRepeats.header

mv TMP/simpleRepeats.bed.gz ./
mv TMP/simpleRepeats.bed.gz.tbi ./


#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">${whatis}</entry>
        <entry key="url"><a>${url}</a></entry>
</properties>
EOF
"""
}

/* repeat masker intervals */
process DOWNLOAD_RMSK {
afterScript "rm -rf TMP"
memory "2g"
input:
	val(meta)
	val(genomeId)
output:
	tuple val("RMSK"),path("rmsk.bed.gz"),path("rmsk.header"),val("CHROM,FROM,TO,RMSK"),val(""),emit:output
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
	def url = genome.rmsk_url
	def whatis="Repeat Masker from ${url}"
"""
hostname 1>&2
mkdir -p TMP
${moduleLoad("htslib jvarkit bedtools")}

set -o pipefail
wget -O - "${url}" |\
	gunzip -c |\
	cut -f6-8 |\
	java -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${reference}" --column 1 --convert SKIP  |\
		sort  -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n |\
		bedtools merge |\
		sed 's/\$/\t1/' |\
		bgzip > TMP/rmsk.bed.gz && \
	tabix -p bed -f TMP/rmsk.bed.gz

mv TMP/rmsk.bed.gz ./
mv TMP/rmsk.bed.gz.tbi ./

echo '##INFO=<ID=RMSK,Number=0,Type=Flag,Description="${whatis}">' > rmsk.header

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">${whatis}</entry>
        <entry key="url"><a>${url}</a></entry>
</properties>
EOF
"""
}

process FREQUENT_GNOMADSV {
afterScript "rm -rf TMP"
input:
	val(meta)
	val(tabix)
output:
	tuple val("GNOMADSV"),path("gnomad.sv.bed.gz"),path("gnomad.sv.header"),val("CHROM,FROM,TO,GNOMAD_SV"),val("GNOMAD_SV:unique"),emit:output
        path("version.xml"),emit:version
script:
	def whatis ="Gnomad SV"
	def pop = params.annotations.gnomad_sv_pop?:"POPMAX_AF"
"""
hostname 1>&2
${moduleLoad("htslib")}
set -o pipefail
mkdir -p TMP

gunzip -c  '${tabix}' |\
	awk '(NR==1) {C=-1;for(i=1;i<=NF;i++) if(\$i=="${pop}") C=i;next;} {if(\$C!="NA" && \$C*1.0 > ${params.annotations.gnomad_sv_AF} ) printf("%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,\$4);}'  |\
	sort -T . -t '\t' -k1,1 -k2,2n |\
	bgzip > TMP/gnomad.sv.bed.gz

tabix -p bed -f TMP/gnomad.sv.bed.gz

mv TMP/gnomad.sv.bed.gz ./
mv TMP/gnomad.sv.bed.gz.tbi ./

echo '##INFO=<ID=GNOMAD_SV,Number=.,Type=String,Description="${whatis} ${pop}/AF greater than ${params.annotations.gnomad_sv_AF}  ">' > gnomad.sv.header

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">${whatis}</entry>
</properties>
EOF
"""
}


process DOWNLOAD_VISTA_ENHANCERS {
afterScript "rm -rf TMP"
tag "${genomeId}"
input:
	val(meta)
	val(genomeId)
output:
	tuple val("VISTA"),path("vistaEnhancers.bed.gz"),path("vistaEnhancers.header"),val("CHROM,FROM,TO,VISTA_ENHANCER"),val("VISTA_ENHANCER:unique"),emit:output
        path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]	
	def url = genome.vista_enhancers_url
	def reference = genome.fasta
	def whatis = "Vista Enhancers from ${url}"
"""
hostname 1>&2
${moduleLoad("htslib jvarkit bedtools")}
set -o pipefail
mkdir -p TMP

wget -O - "${url}" |\
	gunzip -c |\
	cut -f2-5 |\
	java -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${reference}" --column 1 --convert SKIP  |\
	sort -t '\t' -k1,1 -k2,2n -T . |\
	bgzip > TMP/vistaEnhancers.bed.gz


tabix -p bed -f TMP/vistaEnhancers.bed.gz

mv TMP/vistaEnhancers.bed.gz ./
mv TMP/vistaEnhancers.bed.gz.tbi ./

echo '##INFO=<ID=VISTA_ENHANCER,Number=.,Type=String,Description="${whatis}">' > vistaEnhancers.header


#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">${whatis}</entry>
        <entry key="url"><a>${url}</a></entry>
</properties>
EOF
"""
}

process DOWNLOAD_REGFEATURES {
afterScript "rm -rf TMP"
input:
	val(meta)
	val(genomeId)
output:
	tuple val("REGFEAT"),path("regulatory_features.bed.gz"),path("regulatory_features.header"),val("CHROM,FROM,TO,REGULATORY_FEAT"),val("REGULATORY_FEAT:unique"),emit:output
        path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
	def url = genome.ensembl_regulatory_gff_url
	def whatis = "Regulatory Features from ${url}"

"""
hostname 1>&2
${moduleLoad("htslib jvarkit bedtools")}
set -o pipefail
mkdir -p TMP

wget -O - "${url}" |\
	gunzip -c |\
	awk -F '\t' '/^[^#]/ {printf("%s\t%d\t%s\t%s\\n",\$1,int(\$4)-1,\$5,\$3);}' |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP  |\
	LC_ALL=C sort -t '\t' -k1,1 -k2,2n -T . |\
	bgzip > TMP/regulatory_features.bed.gz

tabix -p bed -f TMP/regulatory_features.bed.gz

mv TMP/regulatory_features.bed.gz ./
mv TMP/regulatory_features.bed.gz.tbi ./

echo '##INFO=<ID=REGULATORY_FEAT,Number=.,Type=String,Description="${whatis}">' > regulatory_features.header

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">${whatis}</entry>
        <entry key="url"><a>${url}</a></entry>
</properties>
EOF
"""
}

/**

http://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=regulation&hgta_track=encodeCcreCombined&hgta_table=encodeCcreCombined&hgta_doSchema=describe+table+schema

This track displays the ENCODE Registry of candidate cis-Regulatory Elements (cCREs) in the human genome, a total of 926,535 elements identified and classified by the ENCODE Data Analysis Center according to biochemical signatures. 

*/
process DOWNLOAD_ENCODE_CCRE {
afterScript "rm -rf TMP"
input:
	val(meta)
	val(genomeId)
	path(sortedbed)
output:
	tuple val("CCRE"),path("encode.ccre.bed.gz"),path("encode.ccre.header"),val("CHROM,FROM,TO,ENCODE_CCRE"),val("ENCODE_CCRE:unique"),emit:output
        path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
	def url = genome.encode_ccre_url
	def whatis = "ENCODE Registry of candidate cis-Regulatory Elements (cCREs) in the human genome from ${url}"

"""
hostname 1>&2
${moduleLoad("htslib jvarkit bedtools ucsc")}
set -o pipefail
mkdir -p TMP/CACHE


wget -O TMP/jeter.bb "${url}"

bigBedToBed -udcDir=TMP/CACHE TMP/jeter.bb stdout |\
	awk -F '\t' '{printf("%s\t%s\t%s\t%s|%s|%s\\n",\$1,\$2,\$3,\$4,\$5,\$11)}' |\
	java -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${reference}" --column 1 --convert SKIP  |\
	${sortedbed.name.equals("NO_FILE")?"":" sort -T TMP -t '\t' -k1,1 -k2,2n | bedtools intersect -u -a - -b '${sortedbed}' | "} \
	LC_ALL=C sort -t '\t' -k1,1 -k2,2n -T TMP |\
	bgzip > TMP/encode.ccre.bed.gz


tabix -p bed -f TMP/encode.ccre.bed.gz

mv TMP/encode.ccre.bed.gz ./
mv TMP/encode.ccre.bed.gz.tbi ./

echo '##INFO=<ID=ENCODE_CCRE,Number=.,Type=String,Description="format=NAME|SCORE|LABEL    ${whatis}">' > encode.ccre.header

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">${whatis}</entry>
        <entry key="url"><a>${url}</a></entry>
</properties>
EOF
"""
}


process DOWNLOAD_CLINVAR {
afterScript "rm -rf TMP"
memory "3g"
input:
      	val(meta)
        val(genomeId)
output:
       	tuple val("CLINVAR"),path("clinvar.bcf"),path("clinvar.header"),val("CLNSIG,CLN_ALLELEID"),val(""),emit:vcf
        path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
       	def url = genome.clinvar_vcf_url
        def whatis="ClinVar aggregates information about genomic variation and its relationship to human health."

"""
hostname 1>&2
${moduleLoad("bcftools jvarkit")}
set -o pipefail
mkdir -p TMP

wget -O - "${url}" |\
	gunzip -c |\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar vcfsetdict -R "${genome.fasta}"  -n SKIP |\
	bcftools sort -T TMP/tmp -O b -o TMP/clinvar.bcf

bcftools view --header-only TMP/clinvar.bcf | grep "^##INFO" | cut -d '=' -f3 | cut -d ',' -f1 | grep -v '^CLN' | awk '{printf("INFO/%s\tCLN_%s\\n",\$1,\$1);}' > TMP/rename.tsv
bcftools annotate --rename-annots TMP/rename.tsv -O b -o TMP/jeter.bcf TMP/clinvar.bcf
mv TMP/jeter.bcf TMP/clinvar.bcf

bcftools index --force TMP/clinvar.bcf

mv TMP/clinvar.bcf ./
mv TMP/clinvar.bcf.csi ./

#useless because vcf
touch clinvar.header

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">${whatis}</entry>
        <entry key="url"><a>${url}</a></entry>
        <entry key="versions">${getVersionCmd("bcftools")}</entry>
</properties>
EOF

"""
}



