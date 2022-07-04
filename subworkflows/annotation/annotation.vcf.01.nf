include {jvarkit;isHg19;isHg38;isBlank;hasFeature;moduleLoad;getGnomadGenomePath;getGnomadExomePath} from '../../modules/utils/functions.nf'
include {CONCAT_FILES_01} from '../../modules/utils/concat.files.nf'
include {DOWNLOAD_GTF_01} from '../../modules/gtf/download.gtf.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {DOWNLOAD_GFF3_01} from '../../modules/gff3/download.gff3.01.nf'

boolean isHs37d5(reference) {
	return reference.toString().contains("hs37d5");
	}

String getVcfId(meta,reference) {
	if(!isBlank(meta.vcfid)) return meta.vcfid;
	if(isHs37d5(reference)) return "/LAB-DATA/BiRD/resources/species/human/cng.fr/hs37d5/20210605.hs37d5.ID.sites.vcf.gz";
	return "";
	}



String getAmalgamionVcf(meta) {
	if(!isBlank(meta.joinvcf)) return meta.joinvcf;
	if(meta.reference.contains("hs37d5")) return "/LAB-DATA/BiRD/shares/ITX/u1087/AMALGAMION/20210603.hs37d5.amalgamion.genotyped.bcf";
	return "";
	}

String getAmalgamionSuffix(meta) {
	if(!meta.joinsuffix.isEmpty()) return meta.joinsuffix;
	if(isHg19()) return "_AMALGAMION";
	return "";
	}

boolean isSoftFilter(def meta,String v) {
	if(!meta.containsKey("hardfilters")) return true;
	String[] tokens = meta.hardfilters.split("[, ;]+");
	for(String t: tokens) {
		t= t.trim();
		if(t.equals("*")) return false;
		if(t.equalsIgnoreCase(v)) return false;
		}
	return true;
	}



process complementCaptures {
memory "3g"
input:
	meta
output:
	path("norm.captures.tsv"),emit:output
script:
if(isBlank(meta.captures))
"""
touch norm.captures.tsv
"""
else
"""
hostname 1>&2
${moduleLoad("bedtools jvarkit htslib")}
set -o pipefail

touch norm.captures.tsv

cut -f1,2 "${meta.reference}.fai" |\
	sort -t '\t' -k1,1 -k2,2n -T .  > jeter.genome

grep -v "^#" "${meta.captures}" | while read CN CF
	do

		(gunzip -c "\${CF}" || cat "\${CF}") | grep -v -E '^(browser|track|#)' | cut -f1,2,3 |\
		java -Xmx${task.memory.giga}g -Djava.io.tmpdir=.  -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${meta.reference}" --column 1 --convert SKIP  |\
               		sort -t '\t' -k1,1 -k2,2n -T . |\
			bedtools merge |\
			bedtools complement -i - -g jeter.genome |\
			awk -F '\t' '(int(\$2) < int(\$3)) {printf("%s\t1\\n",\$0);}' |\
			bedtools sort -faidx "${meta.reference}.fai" |\
               		bgzip > "\${CN}.bed.gz"

			tabix -f -p bed "\${CN}.bed.gz"

			echo "\${CN}\t\${PWD}/\${CN}.bed.gz" >> norm.captures.tsv
	done

rm jeter.genome
"""
}

process allBed {
executor "local"
input:
	val(meta)
output:
	path("merged.bed"),emit:merged_bed
script:
"""
set -o pipefail
xargs -a "${meta.beds}" -n 1 cut -f1,2,3| sort -T . -t '\t' -k1,1 -k2,2n > merged.bed
test -s merged.bed
"""
}

/** build snpeff Database from gtf */
process BUILD_SNPEFF {
afterScript "rm -f org.tsv genes.tsv"
memory "10g"
input:
        val(meta)
        val(reference)
        path(gtf)
output:
	path("data"),emit:output
	path("version.xml"),emit:version
script:
	def dbName = file(reference).getSimpleName()
"""
hostname 1>&2
${moduleLoad("snpEff")}
set -o pipefail

mkdir -p "data/${dbName}"
ln -s "${reference}" "data/${dbName}/sequences.fa"

gunzip -c "${gtf}" > data/${dbName}/genes.gtf

# write snpEff contig
cat << EOF > snpEff.config
data.dir = \${PWD}/data/
${dbName}.genome = Human
${dbName}.reference = ${reference}
EOF

# build database
java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar \${SNPEFF_JAR}  build -gtf22 -v "${dbName}" 2> snpeff.errors

rm snpeff.errors

test -s "data/${dbName}/snpEffectPredictor.bin"

rm  data/${dbName}/genes.gtf

cat <<- EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">build custom snpeff from gtf</entry>
        <entry key="gtf">${gtf}</entry>
        <entry key="fasta">${reference}</entry>
        <entry key="output">\${PWD}/snpEff.config</entry>
        <entry key="snpeff">\$(java -jar \${SNPEFF_JAR} -version)</entry>
</properties>
EOF
"""
}

/* gene symbols as bed file sorted on gene symbol */
process BED_FOR_GENES {
afterScript "rm -f org.tsv genes.tsv"
input:
	val(meta)
	val(reference)
      	path(gtf)
output:
       	path("genes.symbols.bed"),emit:bed
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bedtools jvarkit")}
set -o pipefail

java -jar ${jvarkit("gtf2bed")} --columns "gtf.feature,gene_name" -R "${reference}" "${gtf}" |\
awk -F '\t' '(\$4=="gene")' |\
cut -f1,2,3,5 |\
sort -T . -t '\t' -k4,4 > genes.symbols.bed


#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
	<entry key="description">convert gtf to bed</entry>
        <entry key="gtf">${gtf}</entry>
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
        <entry key="url">${url}</entry>
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
        <entry key="url">${url}</entry>
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
        <entry key="url">${url}</entry>
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


wget -O - "${url}" |\
	gunzip -c |\
	grep -v '^chromosome' |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP  |\
	awk -F '\t' '{for(i=1;i<=NF;i++) {printf("%s%s",(i==1?"":(i<5?"\t":"|")),\$i);} printf("\\n");}' |\
	bgzip > greendb.bed.gz

tabix -p bed greendb.bed.gz

echo '##INFO=<ID=GREENDB,Number=.,Type=String,Description="${whatis}. ${url} Format: (regionID|std_type|DB_source|N_sources|N_methods|constrain_pct|PhyloP100_median|closestGene_symbol|closestGene_dist|closestProt_symbol|closestProt_dist|closestGene_TSS_symbol|closestGene_TSS_dist|closestProt_TSS_symbol|closestProt_TSS_dist|controlled_genes|N_controlled_genes)">' > greendb.header

echo "\${PWD}/greendb.bed.gz\tCHROM,FROM,TO,GREENDB\t\${PWD}/greendb.header" > greendb.tsv

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
	<entry key="description">${whatis}</entry>
        <entry key="url">${url}</entry>
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
        <entry key="url">${url}</entry>
</properties>
EOF

"""
else
"""
touch organizer.tsv
echo "<properties/>" > version.xml
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
        <entry key="obo.url">${url0}</entry>
        <entry key="annot.url">${url1}</entry>
</properties>
EOF

"""
else
"""
touch go.terms.tsv
touch go.tsv
echo '<properties/>' > version.xml
"""
}

process DOWNLOAD_TOMMO {
afterScript "rm -f tmp.bed"
input:
	val(meta)
	val(reference)
      	path bed from merged_bed
output:
       	path("tommo.bed.gz") into tommo_vcf
        path("tommo.bed.gz.tbi")
script:
if(isHg19(reference) && hasFeature(meta,"tommo"))
"""
hostname 1>&2
module load htslib/0.0.0 jvarkit bedtools/0.0.0 bcftools/0.0.0

for U in "tommo-14KJPN-20211208-GRCh37_lifted_from_GRCh38-af-autosome" 
do
	wget -O - "https://jmorp.megabank.tohoku.ac.jp/dj1/datasets/tommo-14kjpn-20211208-af_snvindelall/files/\${U}.vcf.gz?download=true" |\
		gunzip -c|\
		java -jar \${JVARKIT_DIST}/vcfsetdict.jar -R "${meta.reference}" |\
		bcftools norm -m '-' --targets-file "${bed}" -O u |\
		bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\n' >> tmp.bed
done
LC_ALL=C sort -t '\t' -k1,1 -k2,2n -T . tmp.bed | uniq |  bgzip > tommo.bed.gz
tabix -s 1 -b 2 -e 2 tommo.bed.gz
"""
else
"""
module load htslib/0.0.0
touch tommo.bed
bgzip tommo.bed
tabix -s 1 -b 2 -e 2 tommo.bed.gz
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
        <entry key="url">${url}</entry>
</properties>
EOF

"""
else
"""
touch remap.tsv
echo "<properties/>" > version.xml
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
        <entry key="url">${url}</entry>
</properties>
EOF


"""
else
"""
touch simpleRepeats.tsv
echo "<properties/>" > version.xml
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
        <entry key="url">${url}</entry>
</properties>
EOF

"""
else
"""
touch rmsk.tsv
echo "<properties/>" > version.xml
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

if(!isBlank(url))
"""
hostname 1>&2
${moduleLoad("htslib jvarkit bedtools")}
set -o pipefail

wget -q  -O - "${url}" |\
	gunzip -c  |\
	awk '(NR==1) {C=-1;for(i=1;i<=NF;i++) if(\$i=="POPMAX_AF") C=i;next;} {if(\$C!="NA" && \$C*1.0 > ${meta.gnomadSVAF} ) printf("%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,\$C);}'  |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${meta.reference}" --column 1 --convert SKIP |\
	sort -T . -t '\t' -k1,1 -k2,2n |\
	bgzip > gnomad.sv.bed.gz

tabix -p bed -f gnomad.sv.bed.gz

echo '##INFO=<ID=GNOMAD_SV,Number=.,Type=String,Description="${whatis}">' > gnomad.sv.header

echo "\${PWD}/gnomad.sv.bed.gz\tCHROM,FROM,TO,GNOMAD_SV\t\${PWD}/gnomad.sv.header" > gnomad.sv.tsv

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">${whatis}</entry>
        <entry key="url">${url}</entry>
</properties>
EOF

"""
else
"""
touch gnomad.sv.tsv
echo '<properties/>' > version.xml
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
	def url = "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/vistaEnhancers.txt.gz"
	def whatis = "Vista Enhancers from ${url}"
if(isHg19(reference) && hasFeature(meta, "vista"))
"""
hostname 1>&2
module load htslib/0.0.0 jvarkit bedtools/0.0.0
set -o pipefail

wget -O - "${url}" |\
	gunzip -c |\
	cut -f2-5 |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${meta.reference}" --column 1 --convert SKIP  |\
	sort -t '\t' -k1,1 -k2,2n -T . |\
	bgzip > vistaEnhancers.bed.gz

tabix -p bed -f vistaEnhancers.bed.gz

echo '##INFO=<ID=VISTA_ENHANCER,Number=.,Type=String,Description="${whatis}">' > vistaEnhancers.header

echo "\${PWD}/vistaEnhancers.bed.gz\tCHROM,FROM,TO,VISTA_ENHANCER\t\${PWD}/vistaEnhancers.header" > vistaEnhancers.tsv

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">${whatis}</entry>
        <entry key="url">${url}</entry>
</properties>
EOF

"""
else
"""
touch "vistaEnhancers.tsv"
echo "<properties/>" > version.xml
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
	def url = "http://ftp.ensembl.org/pub/grch37/release-99/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20180925.gff.gz"
	def whatis = "Regulatory Features from ${url}"

if(isHg19(reference) && hasFeature(meta, "regulatory_features"))
"""
hostname 1>&2
${moduleLoad("htslib jvarkit bedtools")}
set -o pipefail

wget -O - "${url}" |\
	gunzip -c |\
	awk -F '\t' '/^[^#]/ {printf("%s\t%d\t%s\t%s\\n",\$1,int(\$4)-1,\$5,\$3);}' |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${meta.reference}" --column 1 --convert SKIP  |\
	LC_ALL=C sort -t '\t' -k1,1 -k2,2n -T . |\
	bgzip > regulatory_features.bed.gz

tabix -p bed -f regulatory_features.bed.gz

echo '##INFO=<ID=REGULATORY_FEAT,Number=.,Type=String,Description="${whatis}">' > regulatory_features.header

echo "\${PWD}/regulatory_features.bed.gz\tCHROM,FROM,TO,REGULATORY_FEAT\t\${PWD}/regulatory_features.header" > regulatory_features.tsv

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">${whatis}</entry>
        <entry key="url">${url}</entry>
</properties>
EOF
"""
else
"""
touch "regulatory_features.tsv"
echo "<properties/>" > version.xml
"""
}



process APPLY_ANNOTATION {
tag "vcf:${row.vcf} bed:${row.bed}"
cache 'lenient'
afterScript "rm -rf TMP jeter.123.bed jeter1.vcf jeter2.vcf jeter1.bcf jeter1.bcf.csi jeter.tab jeter.tab.gz jeter.tab.gz.tbi hdr.txt"
memory '5 g'
maxForks 30

input:
	val(meta)
	val(reference)
	path(annotations_files)
	path(gtf)
	path(gff3)
	val(row)
output:
	tuple val("${row.bed}"),path("contig.bcf"),emit:bedvcf
	path("contig.bcf.csi"),emit:index
	path("version.xml"),emit:version
script:
	def bed=row.bed
	def vcf=row.vcf
	def cases_file = row.cases?:""
	def controls_file = row.controls?:""
	def dbName = file(reference).getSimpleName()
	def extraBcfTools = meta.extraBcfTools?:""
	def lowGQ = meta.lowGQ?:"70" //TODO a verifier
"""
	hostname 1>&2
	${moduleLoad("bcftools jvarkit snpEff bedtools htslib")}
	set -o pipefail
	set -x
	mkdir -p TMP
	
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Annotation of a VCF file</entry>
		<entry key="bed">${bed}</entry>
		<entry key="vcf">${vcf}</entry>
		<entry key="bcftools.version">\$(bcftools --version | head -n 2 | paste -sd ' ')</entry>
		<entry key="steps">
	EOF


	# normalize bed, bcftools doesn't like more than 3 columns...
	cut -f1,2,3 "${bed}" | bedtools sort -faidx "${reference}.fai" > TMP/jeter.123.bed

	bcftools view ${extraBcfTools} --regions-file "TMP/jeter.123.bed" -O v "${vcf}" |\
		java -jar ${jvarkit("vcfpar")} > TMP/jeter1.vcf


	# samples in pedigree
	if [ ! -z "${!isBlank(cases_file) && !isBlank(controls_file) && hasFeature(meta,"keepSamplesInPed")?"Y":""}" ] ; then
		cat "${cases_file}" "${controls_file}"  | sort -T TMP | uniq > TMP/samples.a
		bcftools query -l TMP/jeter1.vcf  |	sort -T TMP | uniq > TMP/samples.b
		comm -12 TMP/samples.a TMP/samples.b > TMP/samples.c
		test -s TMP/samples.c
		bcftools view --samples-file TMP/samples.c -O u TMP/jeter1.vcf |\
			bcftools view -i 'AC[*]>0' > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi

	# select with jvarkit
	if [ ! -z "${isBlank(meta.select)?"":"Y"}" ] ; then
		java  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/vcffilterjdk.jar --nocode  ${meta.select} TMP/jeter1.vcf > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi
	
	
	# annotate variant ID
	if [ ! -z "${hasFeature(meta,"vcfid")?"Y":""}" ] ; then

		# file must be indexed
		bcftools view -O b -o TMP/jeter1.bcf TMP/jeter1.vcf
		bcftools index TMP/jeter1.bcf
		
		bcftools annotate --annotations "${getVcfId(meta,reference)}"  --regions-file "TMP/jeter.123.bed"  -c ID -O v -o TMP/jeter1.vcf TMP/jeter1.bcf
		rm TMP/jeter1.bcf TMP/jeter1.bcf.csi
	fi


	# polyx
	if [ ! -z "${hasFeature(meta,"polyx") && (meta.polyx as Integer)>0 ?"Y":""}" ] ; then
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar ${jvarkit("vcfpolyx")} --filter ${meta.polyx} --reference ${reference} --tag POLYX TMP/jeter1.vcf > TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi
	

	# allelic ratio
	if [ ! -z "${hasFeature(meta,"AD")?"Y":""}" ] ; then
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar ${jvarkit("vcffilterjdk")} \
		--nocode \
		${isSoftFilter(meta,"AD_RATIO")?"--filter HET_BAD_AD_RATIO ":""} \
		-e 'return variant.getGenotypes().stream().filter(G->G.isHet() && G.hasAD() && G.getAD().length==2).allMatch(G->{int array[]=G.getAD();double r= array[1]/(double)(array[0]+array[1]);return (r>=0.2 && r<=0.8) ;});' TMP/jeter1.vcf > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi

	# GQ
	if [ ! -z "${hasFeature(meta,"GQ")?"Y":""}" ] ; then
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar ${jvarkit("vcffilterjdk")} \
		--nocode \
		${isSoftFilter(meta,"LOW_GQ")?"--filter LOW_GQ${meta.lowGQ} ":""} \
		-e 'return variant.getGenotypes().stream().filter(G->(G.isHet() || G.isHomVar()) && G.hasGQ()).allMatch(G->G.getGQ()>=${meta.lowGQ});' TMP/jeter1.vcf > TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi

	# DP
	if [ ! -z "${hasFeature(meta,"DP")?"Y":""}" ] ; then
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar ${jvarkit("vcffilterjdk")} \
		--nocode \
		${isSoftFilter(meta,"LOW_DP")?"--filter LOW_DP${meta.lowDP} ":""} \
		-e 'return variant.getGenotypes().stream().filter(G->(G.isHet() || G.isHomVar()) && G.hasDP()).allMatch(G->G.getDP()>= ${meta.lowDP} );' TMP/jeter1.vcf > TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi

	# MQ
	if [ ! -z "${hasFeature(meta,"MQ")?"Y":""}" ] ; then
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar ${jvarkit("vcffilterjdk")} \
		--nocode \
		${isSoftFilter(meta,"LOW_MQ")?"--filter LOW_MQ${meta.lowMQ} ":""} \
		-e 'return !variant.hasAttribute("MQ") || variant.getAttributeAsDouble("MQ",1000.0) >= ${meta.lowMQ};' TMP/jeter1.vcf > TMP/jeter2.vcf

	mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi

	# MQRankSum
	if [ ! -z "${hasFeature(meta,"MQRankSum")?"Y":""}" ] ; then
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar ${jvarkit("vcffilterjdk")} \
		--nocode \
		${isSoftFilter(meta,"BAD_MQRankSum")?"--filter BAD_MQRankSum${meta.mqRankSum} ":""} \
		-e 'return !variant.hasAttribute("MQRankSum") || Math.abs(variant.getAttributeAsDouble("MQRankSum",0.0)) <= ${meta.mqRankSum};' TMP/jeter1.vcf > TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi

	# ReadPosRankSum
	if [ ! -z "${hasFeature(meta, "ReadPosRankSum")?"Y":""}" ] ; then
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar ${jvarkit("vcffilterjdk")} \
		--nocode \
		${isSoftFilter(meta,"BAD_ReadPosRankSum")?"--filter BAD_ReadPosRankSum${meta.readPosRankSum} ":""} \
		-e 'return !variant.hasAttribute("ReadPosRankSum") || Math.abs(variant.getAttributeAsDouble("ReadPosRankSum",0.0)) <= ${meta.readPosRankSum};' TMP/jeter1.vcf > TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi

	# SOR
	if [ ! -z "${hasFeature(meta,"SOR")?"Y":""}" ] ; then
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar ${jvarkit("vcffilterjdk")} \
		--nocode \
		${isSoftFilter(meta,"BAD_SOR")?"--filter BAD_SOR${meta.sor} ":""} \
		-e 'return !variant.hasAttribute("SOR") || variant.getAttributeAsDouble("SOR",10000.0) <= ${meta.sor};' TMP/jeter1.vcf > TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi


	# ANNOTATIONS
	cat "${annotations_files}" | while read DATABASE COLS HEADER
	do
		echo "#ANNOT \${DATABASE}" 1>&2
		test -f "\${DATABASE}"
		test -f "\${HEADER}"
		test -f "\${DATABASE}.tbi"

		bcftools annotate --annotations "\${DATABASE}" \
			-h "\${HEADER}" \
			-c "\${COLS}" -O v -o TMP/jeter2.vcf TMP/jeter1.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf

		cat <<- EOF >> version.xml
		<properties>
			<entry key="description">annotate</entry>
			<entry key="header">\${HEADER}</entry>
			<entry key="info.columns">\${COLS}</entry>
		</properties>
		EOF
	done


        if [ ! -z "${gff3.name.equals("NO_FILE")?"":"Y"}" ] ; then

  		bcftools csq -O v --force --local-csq --ncsq 10000 --fasta-ref "${reference}" --gff-annot "${gff3}" TMP/jeter1.vcf > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf

		cat <<- EOF >> version.xml
		<properties>
			<entry key="description">bcftools csq</entry>
			<entry key="gff3">${gff3}</entry>
		</properties>
		EOF

	fi


	if [ ! -z "${isHg19(reference) && hasFeature(meta, "snpeff")?"Y":""}"  ] ; then
		java  -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar "\${SNPEFF_JAR}" eff \
			-config  "\${SNPEFF_CONFIG}" \
			-nodownload -noNextProt -noMotif -noInteraction -noLog -noStats -chr chr -i vcf -o vcf "${isHg19(reference)?"GRCh37.75":"TODO"}" TMP/jeter1.vcf > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi


        if [ ! -z "${isHg19(reference) && hasFeature(meta, "vep")?"Y":""}" ] ; then
		module load ensembl-vep/104.3
		vep --cache  --format vcf --force_overwrite --output_file STDOUT --no_stats --offline  --dir_cache /LAB-DATA/BiRD/resources/apps/vep  --species homo_sapiens --cache_version 104 --assembly GRCh37 --merged --fasta "${reference}" --use_given_ref --vcf < TMP/jeter1.vcf > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf
		module unload ensembl-vep/104.3
	fi

        if [ ! -z "${(hasFeature(meta,"snpeff") || hasFeature(meta, "vep")) && !isBlank(meta.soacn) ?"Y":""}" ] ; then
		java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar ${jvarkit("vcffilterso")} \
			${isSoftFilter(meta,"BAD_SO")?"--filterout  BAD_SO":""} \
			--acn "${meta.soacn}" TMP/jeter1.vcf > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi

	## cadd
	if [ ! -z  "${isHg19(reference) && hasFeature(meta, "cadd")?"Y":""}" ] ; then
        	java  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar ${jvarkit("vcfcadd")} \
			--tabix "/LAB-DATA/BiRD/resources/species/human/krishna.gs.washington.edu/download/CADD/v1.6/whole_genome_SNVs.tsv.gz" TMP/jeter1.vcf > TMP/jeter2.vcf
	      	mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi

	if [ ! -z "${hasFeature(meta,"gnomadGenome")?"Y":""}" ] ; then
		java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar ${jvarkit("vcfgnomad")} \
			--bufferSize 10000 \
			--max-af ${meta.gnomadAF} \
			--prefix "${isSoftFilter(meta,"GNOMAD")?"GNOMAD":""}" \
			--gnomad "${getGnomadGenomePath(meta,reference)}" --fields "${meta.gnomadPop}" TMP/jeter1.vcf > TMP/jeter2.vcf
                mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi

	if [ ! -z "${hasFeature(meta, "gnomadExome")?"Y":""}" ] ; then
		java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar ${jvarkit("vcfgnomad")} \
			--bufferSize 10000 \
			--max-af ${meta.gnomadAF} \
			--prefix "${isSoftFilter(meta,"GNOMAD")?"GNOMAD":""}" \
			--gnomad "${getGnomadExomePath(meta,reference)}" --fields "${meta.gnomadPop}" TMP/jeter1.vcf > TMP/jeter2.vcf
                mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi
	



	bcftools sort --max-mem "${task.memory.giga}g" -T TMP -O b -o TMP/contig.bcf TMP/jeter1.vcf
	bcftools index TMP/contig.bcf

	mv TMP/contig.bcf ./
	mv TMP/contig.bcf.csi ./

	cat <<- EOF >> version.xml
	</entry>
	</properties>
	EOF


"""
stub:
"""
	# tommo VCF japan variants
	if [ "${isHg19(refrence) && hasFeature(meta,"tommo")?"Y":"N"}" == "Y" ] ; then

		echo -e '##INFO=<ID=AC_TOMMO,Number=A,Type=Integer,Description="AC from Tommo DB">' > hdr.txt
		echo -e '##INFO=<ID=AF_TOMMO,Number=A,Type=Float,Description="AF from Tommo DB">' >> hdr.txt
		echo -e '##INFO=<ID=AN_TOMMO,Number=1,Type=Integer,Description="AF from Tommo DB">' >> hdr.txt
		
		bcftools annotate  --mark-sites +IN_TOMMO -a "${tommo}" -h hdr.txt -c 'CHROM,POS,REF,ALT,AC_TOMMO,AN_TOMMO,AF_TOMMO' jeter1.vcf > jeter2.vcf
		mv jeter2.vcf jeter1.vcf

		rm hdr.txt
	fi


	# add external info AF,AC,AN e.g. amalgamion
	if [ ! -z "${hasFeature(meta,"joinvcf")?"Y":""}" ] ; then
		# get interval for this vcf
		bcftools query -f '%CHROM\t%POS0\t%END\\n' jeter1.vcf |\
			sort -T . -t '\t' -k1,1 -k2,2n |\
			bedtools merge > jeter.bed

		# if bed is empty
		if ! [ -s  jeter.123.bed ] ; then
			head -n1 "jeter.123.bed" | awk -F '\t' '{printf("%s\t0\t1\\n",\$1);}'  > jeter.bed
		fi

		bcftools norm --regions-file jeter.bed -O u -m- "${getAmalgamionVcf(reference)}" |\
			bcftools  +fill-tags -O u  -- -t AN,AC,AF |\
			bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\\n' > jeter.tab
		rm jeter.bed
		bgzip jeter.tab
		tabix -s1 -b2 -e2 jeter.tab.gz

		echo -e '##INFO=<ID=AC${getAmalgamionSuffix(meta)},Number=A,Type=Integer,Description="AC from ${getAmalgamionVcf()}">' > hdr.txt
		echo -e '##INFO=<ID=AF${getAmalgamionSuffix(meta)},Number=A,Type=Float,Description="AF from ${getAmalgamionVcf()}">' >> hdr.txt
		echo -e '##INFO=<ID=AN${getAmalgamionSuffix(meta)},Number=1,Type=Integer,Description="AF from ${getAmalgamionVcf()}">' >> hdr.txt
		
		bcftools annotate  --mark-sites +IN${getAmalgamionSuffix(meta)} -a jeter.tab.gz -h hdr.txt -c 'CHROM,POS,REF,ALT,AC${getAmalgamionSuffix(meta)},AN${getAmalgamionSuffix(meta)},AF${getAmalgamionSuffix(meta)}' jeter1.vcf > jeter2.vcf
		mv jeter2.vcf jeter1.vcf

		rm hdr.txt jeter.tab.gz jeter.tab.gz.tbi
	fi

	# merge external vcf e.g. amalgamion
	if [ ! -z "${isBlank(mergevcf)?"":"Y"}" ] && [ ! -z "${hasFeature(meta,"amalgamion")?"Y":""}" ] ; then
		# get interval for this vcf
		bcftools query -f '%CHROM\t%POS0\t%END\\n' jeter1.vcf |\
			sort -T . -t '\t' -k1,1 -k2,2n |\
			bedtools merge > jeter.bed

		# if bed is empty
		if ! [ -s  jeter.bed ] ; then
			head -n1 "jeter.123.bed" | awk -F '\t' '{printf("%s\t0\t1\\n",\$1);}'  > jeter.bed
		fi

		# exclude common samples
		bcftools query -l jeter1.vcf |sort | uniq > samples1.txt
		bcftools query -l "${mergevcf}" |sort | uniq > samples2.txt
		comm -13 samples1.txt samples2.txt > samples3.txt
		test -s samples3.txt
		
		if [ -s  samples1.txt ] ; then
			awk '{printf(" %s variant.getGenotype(\\"%s\\").isNoCall()", (NR==1?" return !(" :" && "),\$0);} END {printf(") ;\\n");} ' samples1.txt > jeter.code
		else
			echo "return true;" > jeter.code
		fi

		rm samples1.txt samples2.txt 

		# get variants to join in that region
		bcftools annotate -x 'ID,QUAL,FILTER,INFO' --regions-file jeter.bed -O u "${mergevcf}"  |\
			bcftools view --samples-file samples3.txt -O b -o jeter4.bcf
		rm samples3.txt jeter.bed
		bcftools index jeter4.bcf
	
		# bcftools merge wants indexed file
		bcftools view -O b -o jeter1.vcf.gz jeter1.vcf
		bcftools index jeter1.vcf.gz

		# run merge
		bcftools merge  --regions-file "jeter.123.bed" -O v -o jeter1.vcf jeter1.vcf.gz jeter4.bcf

		# keep where a genotype of original VCF has one called genotype
		java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar ${jvarkit("vcffilterjdk")} -f jeter.code jeter1.vcf > jeter2.vcf
		mv jeter2.vcf jeter1.vcf
		

		bcftools query -f '%CHROM:%POS\\n' jeter1.vcf  | sort | uniq > tmp1.txt
		bcftools query -f '%CHROM:%POS\\n' jeter1.vcf.gz  | sort | uniq > tmp2.txt
		comm -3 tmp1.txt tmp2.txt >&2
		rm tmp1.txt tmp2.txt
		rm jeter4.bcf jeter4.bcf.csi jeter1.vcf.gz jeter1.vcf.gz.csi
		

	fi

	# de novo
	if [ ! -z "${meta.pedigree}" ] ; then
		# check parents
		awk -F '\t' '(\$3!="0" || \$4!="0")' "${meta.pedigree}" > parents.txt
		
		if [ -s parents.txt ] ; then
			java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar ${jvarkit("vcftrio")} \
				--pedigree "${meta.pedigree}" jeter1.vcf > jeter2.vcf
			mv jeter2.vcf jeter1.vcf
		fi
		rm parents.txt

		if [ ! -z "${hasFeature(meta, "contrast")?"Y":""}" ] ; then
			# check cases and controls
			awk -F '\t' '(\$6=="case" || \$6=="1")' "${meta.pedigree}" | cut -f 2 > cases.list
			awk -F '\t' '(\$6=="control" || \$6=="0")' "${meta.pedigree}" | cut -f 2 > ctrl.list

			if [ ! -s TMP/cases.list ] && [ ! -s TMP/ctrl.list ] ; then

    				bcftools +contrast -0 ctrl.list -1 cases.list \
      					-a PASSOC,FASSOC,NASSOC,NOVELAL,NOVELGT -O b -o jeter2.vcf jeter1.vcf
				mv jeter2.vcf jeter1.vcf
			fi
			rm cases.list ctrl.list
		fi
	fi



	
	bcftools sort --max-mem "${task.memory.giga}g" -T . -O b -o contig.bcf jeter1.vcf
	bcftools index contig.bcf
	rm jeter1.vcf
"""
}

workflow ANNOTATE {
	take:
		meta
		reference
		input
	main:
		annot_ch = Channel.empty()
		version_ch = Channel.empty()

		gtf_ch = DOWNLOAD_GTF_01(meta.plus("with_tabix":true,"gtfurl":(isHg19(reference)?"http://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr_patch_hapl_scaff.gtf.gz":"")),reference)
                version_ch = version_ch.mix(gtf_ch.version)

		gff3_ch = DOWNLOAD_GFF3_01(meta.plus("with_tabix":true,"gff3url":(isHg19(reference)?"http://ftp.ensembl.org/pub/grch37/current/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.chr_patch_hapl_scaff.gff3.gz":"")),reference)
		version_ch = version_ch.mix(gff3_ch.version)


		//snpeff_db = BUILD_SNPEFF(meta,reference, gtf_ch.gtf)
                //version_ch = version_ch.mix(snpeff_db.version)

		bed4genes_ch = BED_FOR_GENES(meta, reference, gtf_ch.gtf)
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

		db11_ch = DOWNLOAD_VISTA_ENHANCERS(meta, reference)
		version_ch = version_ch.mix(db11_ch.version)
		annot_ch = annot_ch.mix(db11_ch.output)

		db12_ch = DOWNLOAD_REGFEATURES(meta, reference)
		version_ch = version_ch.mix(db12_ch.version)
		annot_ch = annot_ch.mix(db12_ch.output)
		

		db2_ch = CONCAT_FILES_01(meta,annot_ch.collect())


		annot_vcf_ch = APPLY_ANNOTATION(
			meta,
			reference,
			db2_ch.output,
			gtf_ch.gtf,
			gff3_ch.gff3,
			input
			)
		version_ch = version_ch.mix(annot_vcf_ch.version.first())

		to_file_ch = COLLECT_TO_FILE_01([:],annot_vcf_ch.bedvcf.map{T->T[0]+"\t"+T[1].toRealPath()}.collect())
		version_ch = version_ch.mix(to_file_ch.version)

		version_ch = MERGE_VERSION(meta, "annotation", "VCF annotation", version_ch.collect())
	emit:
		version= version_ch
		bedvcf = to_file_ch.output
	}
