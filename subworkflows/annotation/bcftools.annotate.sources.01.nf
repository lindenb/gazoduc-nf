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

include {getVersionCmd;jvarkit;isHg19;isHg38;isBlank;hasFeature;moduleLoad;getGnomadGenomePath;getGnomadExomePath} from '../../modules/utils/functions.nf'
include {CONCAT_FILES_01} from '../../modules/utils/concat.files.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'




workflow BCFTOOLS_ANNOTATE_SOURCES {
	take:
		meta
		reference
		gtf /** gtf */
		bed /** limit to that bed or NO_FILE */
	main:
		version_ch = Channel.empty()
		annot_ch  = Channel.empty()

		if(bed.name.equals("NO_FILE")) {
			sorted_bed = bed
			}
		else
			{
			sortbed_ch = SORT_BED(meta, reference, bed)
			version_ch = version_ch.mix(sortbed_ch.version)
			sorted_bed = sortbed_ch.bed
			}

		bed4genes_ch = BED_FOR_GENES(meta, reference, gtf_ch.gtf, sorted_bed)
                version_ch = version_ch.mix(bed4genes_ch.version)

		if(hasFeature(meta,"greendb")) {
			db1_ch = DOWNLOAD_GREENDB(meta,reference)
			version_ch = version_ch.mix(db1_ch.version)
			annot_ch = annot_ch.mix(db1_ch.output)
		}

		if(hasFeature(meta,"bhfucl")) {
			db2_ch = DOWNLOAD_BHFUCL(meta, reference, bed4genes_ch.bed)
			version_ch = version_ch.mix(db1_ch.version)
			annot_ch = annot_ch.mix(db1_ch.output)
		}

		if(hasFeature(meta,"diseases")) {
			db3_ch = DOWNLOAD_DISEASES(meta, reference, bed4genes_ch.bed)
			version_ch = version_ch.mix(db3_ch.version)
			annot_ch = annot_ch.mix(db3_ch.output)
		}

		if(hasFeature(meta,"gencc")) {
			db4_ch = DOWNLOAD_GENCC(meta, reference, bed4genes_ch.bed)
			version_ch = version_ch.mix(db4_ch.version)
			annot_ch = annot_ch.mix(db4_ch.output)
		}

		if(hasFeature(meta,"organizer")) {
			db5_ch = DOWNLOAD_ORGANIZER(meta, reference, bed4genes_ch.bed)
			version_ch = version_ch.mix(db5_ch.version)
			annot_ch = annot_ch.mix(db5_ch.output)
		}

		if(hasFeature(meta,"go")) {
			db6_ch = DOWNLOAD_GO(meta,reference, gff3_ch.gff3)
			version_ch = version_ch.mix(db6_ch.version)
			annot_ch = annot_ch.mix(db6_ch.output)
		}

		if(hasFeature(meta, "remap")) {
			db7_ch = DOWNLOAD_REMAP(meta, reference)
			version_ch = version_ch.mix(db7_ch.version)
			annot_ch = annot_ch.mix(db7_ch.output)
		}

		if(hasFeature(meta, "simple_repeats")) {
			db8_ch = DOWNLOAD_SIMPLE_REPEATS(meta, reference)
			version_ch = version_ch.mix(db8_ch.version)
			annot_ch = annot_ch.mix(db8_ch.output)
		}

		if(hasFeature(meta, "rmsk")) {
			db9_ch = DOWNLOAD_RMSK(meta, reference)
			version_ch = version_ch.mix(db9_ch.version)
			annot_ch = annot_ch.mix(db9_ch.output)
		}

		if(hasFeature(meta, "gnomad_sv")) {
			db10_ch = DOWNLOAD_GNOMADSV(meta, reference)
			version_ch = version_ch.mix(db10_ch.version)
			annot_ch = annot_ch.mix(db10_ch.output)
		}

		if(hasFeature(meta, "vista")) {
			db11_ch = DOWNLOAD_VISTA_ENHANCERS(meta, reference)
			version_ch = version_ch.mix(db11_ch.version)
			annot_ch = annot_ch.mix(db11_ch.output)
		}

		if(hasFeature(meta, "regfeatures")) {
			db12_ch = DOWNLOAD_REGFEATURES(meta, reference)
			version_ch = version_ch.mix(db12_ch.version)
			annot_ch = annot_ch.mix(db12_ch.output)
		}


		if(hasFeature(meta, "encode_ccre")) {
			db13_ch = DOWNLOAD_ENCODE_CCRE(meta, reference, sorted_bed)
			version_ch = version_ch.mix(db13_ch.version)
			annot_ch = annot_ch.mix(db13_ch.output)
		}

		db2_ch = CONCAT_FILES_01(meta,annot_ch.collect())

		version_ch = MERGE_VERSION(meta, "bcftoolsAnnotation", "annotations for bcftools", version_ch.collect())
		
	emit:
		version = version_ch
		output = db2_ch.output
	}


process SORT_BED {
tag "${bed}"
input:
	val(meta)
	val(reference)
	path(bed)
output:
	path("sorted.bed"),emit:bed
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bedtools jvarkit")}
set -o pipefail

${bed.name.endsWith(".gz")?"gunzip -c ":"cat "} "${bed}" |\
	grep -vE "^(#|browser|track)" |\
	cut -f1,2,3 |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP  |\
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
	val(reference)
      	path(gtf)
	path(sortedbed)
output:
       	path("genes.symbols.bed"),emit:bed
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bedtools jvarkit")}
set -o pipefail
mkdir -p TMP


java -jar \${JVARKIT_DIST}/gtf2bed.jar --columns "gtf.feature,gene_name" -R "${reference}" "${gtf}" |\
awk -F '\t' '(\$4=="gene")' |\
cut -f1,2,3,5 |\
${sortedbed.name.equals("NO_FILE")?"":" sort -T TMP -t '\t' -k1,1 -k2,2n | bedtools intersect -u -a - -b '${sortedbed}' | "} \
sort -T TMP -t '\t' -k4,4 | uniq > TMP/genes.symbols.bed

mv TMP/genes.symbols.bed ./

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
afterScript "rm -f genes.txt"
input:
	val(meta)
	val(reference)
      	path(bed)
output:
	path("bhfucl.tsv"),emit:output
        path("version.xml"),emit:version
when:
script:
	def url =  "http://ftp.ebi.ac.uk/pub/databases/GO/goa/bhf-ucl/gene_association.goa_bhf-ucl.gz"
	def whatis= "Cardiovascular Gene Ontology Annotation Initiative https://www.ebi.ac.uk/GOA/CVI"

"""
hostname 1>&2
${moduleLoad("htslib")}
set -o pipefail


wget -O - "${url}" |\
	gunzip -c |\
	awk -F '\t' '/#/{next} (\$4 ~ /^NOT/) {next;} {print \$3;}' |\
	uniq | sort | uniq |
	sed 's/\$/\t1/' > genes.txt

join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,2.1,2.2' "${bed}" genes.txt |\
sort -T . -t '\t' -k1,1 -k2,2n | bgzip > bhfucl.bed.gz

tabix -p bed bhfucl.bed.gz


echo '##INFO=<ID=BHFUCL,Number=0,Type=Flag,Description="${whatis}">' > bhfucl.header

echo "\${PWD}/bhfucl.bed.gz\tCHROM,FROM,TO,BHFUCL\t\${PWD}/bhfucl.header" > bhfucl.tsv

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
afterScript "rm -f org.tsv"
input:
	val(meta)
	val(reference)
      	path(bed)
output:
	path("diseases.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def url = "https://download.jensenlab.org/human_disease_knowledge_filtered.tsv"
	def whatis="DISEASES is a weekly updated web resource that integrates evidence on disease-gene associations from automatic text mining, manually curated literature, cancer mutation data, and genome-wide association studies. https://diseases.jensenlab.org/"

if(!isBlank(url))
"""
hostname 1>&2
${moduleLoad("htslib")}
set -o pipefail


wget -O - "https://download.jensenlab.org/human_disease_knowledge_filtered.tsv" |\
	cut -f 2,4 | tr -c "[a-zA-Z_0-9\\-\\n\\t.]" "_" | tr -s "_"  |\
	sort -T . -t '\t' -k1,1  | uniq > org.tsv

join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,2.2' "${bed}" org.tsv |\
sort -T . -t '\t' -k1,1 -k2,2n | bgzip > diseases.bed.gz

tabix -p bed diseases.bed.gz

echo '##INFO=<ID=DISEASES,Number=.,Type=String,Description="${whatis}">' > diseases.header

echo "\${PWD}/diseases.bed.gz\tCHROM,FROM,TO,DISEASES\t\${PWD}/diseases.header" > diseases.tsv

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
input:
	val(meta)
	val(reference)
      	path(bed)
output:
       	path("gencc.tsv"),emit:output
        path("version.xml"),emit:version
script:
	def url= "https://search.thegencc.org/download/action/submissions-export-tsv"
	def whatis = "gencc The Gene Curation Coalition brings together groups engaged in the evaluation of gene-disease validity with a willingness to share data publicly https://thegencc.org/about.html ${url} "

if(!isBlank(url))
"""
hostname 1>&2
${moduleLoad("htslib bedtools")}
set -o pipefail

wget --no-check-certificate -q  -O - "${url}" |\
	awk -F '\t' '(NR==1) {C=-1;D=-1;for(i=1;i<=NF;i++) {if(\$i=="\\"gene_symbol\\"") C=i;if(\$i=="\\"disease_title\\"") D=i;}next;} {if(C<1 ||  D<1) next; G=\$C;H=\$D; gsub(/[^A-Za-z0-9\\.\\-]+/,"_",G);gsub(/[^A-Za-z0-9\\.\\-]+/,"_",H);  printf("%s\t%s\\n",G,H);}'  |\
	sort -T . -t '\t' -k1,1  | uniq > org.tsv

join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,2.2' "${bed}" org.tsv |\
sort -T . -t '\t' -k1,1 -k2,2n | bgzip > gencc.bed.gz

tabix -p bed -f gencc.bed.gz
rm org.tsv

echo '##INFO=<ID=GENCC,Number=.,Type=String,Description="${whatis}">' > gencc.header

echo "\${PWD}/gencc.bed.gz\tCHROM,FROM,TO,GENCC\t\${PWD}/gencc.header" > gencc.tsv

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
afterScript "rm -f org.tsv"
input:
	val(meta)
	val(reference)
	path(sortedbed)
output:
	path("greendb.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def url = (isHg19(reference)?"https://zenodo.org/record/5636209/files/GRCh37_GREEN-DB.bed.gz?download=1":
		   (isHg38(reference)?"https://zenodo.org/record/5636209/files/GRCh38_GREEN-DB.bed.gz?download=1":""))
	def whatis = "Data from greendb https://zenodo.org/record/5636209.  A comprehensive collection of 2.4 million regulatory elements in the human genome collected from previously published databases, high-throughput screenings and functional studies."

if(!isBlank(url))
"""
hostname 1>&2
${moduleLoad("htslib bedtools jvarkit")}
set -o pipefail

if ${!isBlank(url)} ; then

wget -O - "${url}" |\
	gunzip -c |\
	grep -v '^chromosome' |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP  |\
	awk -F '\t' '{for(i=1;i<=NF;i++) {printf("%s%s",(i==1?"":(i<5?"\t":"|")),\$i);} printf("\\n");}' |\
	${sortedbed.name.equals("NO_FILE")?"":" sort -T . -t '\t' -k1,1 -k2,2n | bedtools intersect -u -a - -b '${sortedbed}' | "} \
	sort -T . -t '\t' -k1,1 -k2,2n |\
	bgzip > greendb.bed.gz

else
	touch greendb.bed
	bgzip greendb.bed

fi

tabix -p bed greendb.bed.gz

echo '##INFO=<ID=GREENDB,Number=.,Type=String,Description="${whatis}. ${url} Format: (regionID|std_type|DB_source|N_sources|N_methods|constrain_pct|PhyloP100_median|closestGene_symbol|closestGene_dist|closestProt_symbol|closestProt_dist|closestGene_TSS_symbol|closestGene_TSS_dist|closestProt_TSS_symbol|closestProt_TSS_dist|controlled_genes|N_controlled_genes)">' > greendb.header

echo "\${PWD}/greendb.bed.gz\tCHROM,FROM,TO,GREENDB\t\${PWD}/greendb.header" > greendb.tsv

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
	<entry key="description">${whatis}</entry>
        <entry key="url"><a>${url}</a></entry>
        <entry key="versions">${getVersionCmd("tabix wget jvarkit/bedrenamechr awk")}</entry>
</properties>
EOF
"""
else
"""
touch greendb.tsv
echo '<properties/>' > version.xml
"""
}


process DOWNLOAD_ORGANIZER {
afterScript "rm -f org.tsv"
input:
	val(meta)
	val(reference)
      	path(bed)
output:
	path("organizer.tsv"),emit:output
	path("version.xml"),emit:version
script:

	def url = "http://geneorganizer.huji.ac.il/assets/Uploads/diseases-organs.v12-1.csv.gz"
	def whatis = "Location from Gene Organizer ${url}"

if(!isBlank(url))
"""
hostname 1>&2
${moduleLoad("htslib")}
set -o pipefail


wget -O - "${url}" | gunzip -c |\
cut -d, -f1,3 | uniq |\
tr -d '"' | tr "," "\t" |\
tr " " "_" | sort -t '\t' -T . -k1,1 > org.tsv

join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,2.2' "${bed}" org.tsv |\
sort -T . -t '\t' -k1,1 -k2,2n | bgzip > organizer.bed.gz

tabix -p bed organizer.bed.gz


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
afterScript "rm -f go-basic.obo goa_human.gaf.gz jeter.gff jeter.terms"
memory "10g"
input:
	val(meta)
	val(reference)
      	path(gff3)
output:
       	path("go.tsv"),emit:output
	path("go.terms.tsv"),emit:goterms_tsv
        path("version.xml"),emit:version
script:
	def url0 = "http://purl.obolibrary.org/obo/go/go-basic.obo"
	def url1 = "http://geneontology.org/gene-associations/goa_human.gaf.gz"
	def whatis = "GeneOntology for terms ${meta.goTerms}"

if(!isBlank(meta.goTerms))
"""
hostname 1>&2
set -o pipefail
${moduleLoad("jvarkit bedtools htslib")}

wget -O "go-basic.obo" "${url0}"
wget -O "goa_human.gaf.gz" "${url1}"

echo "${meta.goTerms}" | tr " ,|;" "\\n" | sort |uniq | grep -v '^\$' > jeter.terms
test -s jeter.terms


java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=. -jar ${jvarkit("goutils")} \
	--action dump_table \
	--accession-file "jeter.terms" \
	-go go-basic.obo  | grep -F -f jeter.terms > go.terms.tsv

gunzip -c "${gff3}" | awk -F '\t' '(\$1 ~ /^#/ || \$3=="gene")' | \
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=. -jar ${jvarkit("goutils")} \
	--action gff3 \
	--accession-file "jeter.terms" \
	-go go-basic.obo \
	-goa goa_human.gaf.gz  > jeter.gff

awk '/^#/ {next} {printf("%s\t%d\t%s\\n",\$1,int(\$4)-1,\$5);}' jeter.gff |\
	sort -T . -t '\t' -k1,1 -k2,2n |\
	bedtools merge |\
	sed 's/\$/\t1/' |\
	bgzip > go.bed.gz

tabix -p bed go.bed.gz


echo '##INFO=<ID=GO,Number=0,Type=Flag,Description="${whatis}">' > go.header

echo "\${PWD}/go.bed.gz\tCHROM,FROM,TO,GO\t\${PWD}/go.header" > go.tsv

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
input:
	val(meta)
	val(reference)
output:
	path("remap.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def url = isHg19(reference)?"https://remap.univ-amu.fr/storage/remap2022/hg19/MACS2/remap2022_crm_macs2_hg19_v1_0.bed.gz":""
	def whatis="ReMap is a large scale integrative analysis of DNA-binding experiments for Homo sapiens, Mus musculus, Drosophila melanogaster and Arabidopsis thaliana transcriptional regulators. The catalogues are the results of the manual curation of ChIP-seq, ChIP-exo, DAP-seq from public sources (GEO, ENCODE, ENA). ${url}"

if(!isBlank(url))
"""
hostname 1>&2
${moduleLoad("htslib jvarkit bedtools")}

set -o pipefail
wget -O - "${url}" |\
	gunzip -c |\
	cut -f1,2,3 |\
	sed 's/\$/\t1/' |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${meta.reference}" --column 1 --convert SKIP  |\
		sort -t '\t' -k1,1 -k2,2n -T . |\
		uniq |\
		bgzip > remap.bed.gz && \
	tabix -p bed -f remap.bed.gz


echo '##INFO=<ID=REMAP,Number=.,Type=String,Description="${whatis}">' > remap.header

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
input:
	val(meta)
	val(reference)
output:
	path("simpleRepeats.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def url = isHg19(reference)?"https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz":""
	def whatis = "Simple repeats from ${url}"

if(!isBlank(url))
"""
hostname 1>&2
${moduleLoad("htslib jvarkit bedtools")}
set -o pipefail
wget -O - "${url}" |\
	gunzip -c |\
	cut -f2-4 |\
	sed 's/\$/\t1/' |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${meta.reference}" --column 1 --convert SKIP  |\
		sort -t '\t' -k1,1 -k2,2n -T . |\
		uniq |\
		bgzip > simpleRepeats.bed.gz && \
	tabix -p bed -f simpleRepeats.bed.gz


echo '##INFO=<ID=REPEAT,Number=0,Type=Flag,Description="${whatis}">' > simpleRepeats.header

echo "\${PWD}/simpleRepeats.bed.gz\tCHROM,FROM,TO,REPEAT\t\${PWD}/simpleRepeats.header" > simpleRepeats.tsv

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
memory "2g"
input:
	val(meta)
	val(reference)
output:
	path("rmsk.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def url = isHg19(reference)?"https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz" :""
	def whatis="Repeat Masker from ${url}"

if(!isBlank(url))
"""
hostname 1>&2

${moduleLoad("htslib jvarkit bedtools")}

set -o pipefail
wget -O - "${url}" |\
	gunzip -c |\
	cut -f6-8 |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${meta.reference}" --column 1 --convert SKIP  |\
		sort  -S ${task.memory.kilo} -T . -t '\t' -k1,1 -k2,2n |\
		bedtools merge |\
		sed 's/\$/\t1/' |\
		bgzip > rmsk.bed.gz && \
	tabix -p bed -f rmsk.bed.gz


echo '##INFO=<ID=RMSK,Number=0,Type=Flag,Description="${whatis}">' > rmsk.header

echo "\${PWD}/rmsk.bed.gz\tCHROM,FROM,TO,RMSK\t\${PWD}/rmsk.header" > rmsk.tsv

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

process DOWNLOAD_GNOMADSV {
input:
	val(meta)
	val(reference)
output:
       	path("gnomad.sv.tsv"),emit:output
        path("version.xml"),emit:version
script:
	def url=isHg19(reference)?"https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.bed.gz":""
	def whatis ="Gnomad SV from ${url}"

"""
hostname 1>&2
${moduleLoad("htslib jvarkit bedtools")}
set -o pipefail

if ${!isBlank(url)} ; then

wget -q  -O - "${url}" |\
	gunzip -c  |\
	awk '(NR==1) {C=-1;for(i=1;i<=NF;i++) if(\$i=="POPMAX_AF") C=i;next;} {if(\$C!="NA" && \$C*1.0 > ${meta.gnomadSVAF} ) printf("%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,\$C);}'  |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${meta.reference}" --column 1 --convert SKIP |\
	sort -T . -t '\t' -k1,1 -k2,2n |\
	bgzip > gnomad.sv.bed.gz

else
	touch gnomad.sv.bed
	bgzip gnomad.sv.bed
fi

tabix -p bed -f gnomad.sv.bed.gz

echo '##INFO=<ID=GNOMAD_SV,Number=.,Type=String,Description="${whatis}">' > gnomad.sv.header

echo "\${PWD}/gnomad.sv.bed.gz\tCHROM,FROM,TO,GNOMAD_SV\t\${PWD}/gnomad.sv.header" > gnomad.sv.tsv

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


process DOWNLOAD_VISTA_ENHANCERS {
input:
	val(meta)
	val(reference)
output:
       	path("vistaEnhancers.tsv"),emit:output
        path("version.xml"),emit:version
script:
	def url = isHg19(reference)?"https://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/vistaEnhancers.txt.gz":""
	def whatis = "Vista Enhancers from ${url}"
"""
hostname 1>&2
${moduleLoad("htslib jvarkit bedtools")
set -o pipefail


if ${!isBlank(url)} ; then

wget -O - "${url}" |\
	gunzip -c |\
	cut -f2-5 |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${meta.reference}" --column 1 --convert SKIP  |\
	sort -t '\t' -k1,1 -k2,2n -T . |\
	bgzip > vistaEnhancers.bed.gz

else
	touch vistaEnhancers.bed
	bgzip vistaEnhancers.bed

fi

tabix -p bed -f vistaEnhancers.bed.gz

echo '##INFO=<ID=VISTA_ENHANCER,Number=.,Type=String,Description="${whatis}">' > vistaEnhancers.header

echo "\${PWD}/vistaEnhancers.bed.gz\tCHROM,FROM,TO,VISTA_ENHANCER\t\${PWD}/vistaEnhancers.header" > vistaEnhancers.tsv

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
input:
	val(meta)
	val(reference)
output:
       	path("regulatory_features.tsv"),emit:output
        path("version.xml"),emit:version
script:
	def url = isHg19(reference)?"http://ftp.ensembl.org/pub/grch37/release-99/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20180925.gff.gz":""
	def whatis = "Regulatory Features from ${url}"

"""
hostname 1>&2
${moduleLoad("htslib jvarkit bedtools")}
set -o pipefail

if ${!isBlank(url)} ; then

wget -O - "${url}" |\
	gunzip -c |\
	awk -F '\t' '/^[^#]/ {printf("%s\t%d\t%s\t%s\\n",\$1,int(\$4)-1,\$5,\$3);}' |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP  |\
	LC_ALL=C sort -t '\t' -k1,1 -k2,2n -T . |\
	bgzip > regulatory_features.bed.gz

else

	touch regulatory_features.bed
	bgzip regulatory_features.bed
fi


tabix -p bed -f regulatory_features.bed.gz

echo '##INFO=<ID=REGULATORY_FEAT,Number=.,Type=String,Description="${whatis}">' > regulatory_features.header

echo "\${PWD}/regulatory_features.bed.gz\tCHROM,FROM,TO,REGULATORY_FEAT\t\${PWD}/regulatory_features.header" > regulatory_features.tsv

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
	val(reference)
	path(sortedbed)
output:
       	path("encode.ccre.tsv"),emit:output
        path("version.xml"),emit:version
script:
	def url = isHg38(reference)?"http://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb":""
	def whatis = "ENCODE Registry of candidate cis-Regulatory Elements (cCREs) in the human genome from ${url}"

"""
hostname 1>&2
${moduleLoad("htslib jvarkit bedtools ucsc")}
set -o pipefail
mkdir -p TMP

if ${!isBlank(url)} ; then

wget -O TMP/jeter.bb "${url}"

bigBedToBed TMP/jeter.bb stdout |\
	awk -F '\t' '{printf("%s\t%s\t%s\t%s|%s|%s\\n",\$1,\$2,\$3,\$4,\$5,\$11)}' |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP  |\
	${sortedbed.name.equals("NO_FILE")?"":" sort -T TMP -t '\t' -k1,1 -k2,2n | bedtools intersect -u -a - -b '${sortedbed}' | "} \
	LC_ALL=C sort -t '\t' -k1,1 -k2,2n -T TMP |\
	bgzip > encode.ccre.bed.gz

else

	touch encode.ccre.bed
	bgzip encode.ccre.bed
fi


tabix -p bed -f encode.ccre.bed.gz

echo '##INFO=<ID=ENCODE_CCRE,Number=.,Type=String,Description="format=NAME|SCORE|LABEL    ${whatis}">' > encode.ccre.header

echo "\${PWD}/encode.ccre.bed.gz\tCHROM,FROM,TO,ENCODE_CCRE\t\${PWD}/encode.ccre.header" > encode.ccre.tsv

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
