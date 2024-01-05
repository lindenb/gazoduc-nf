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

include {getVersionCmd;jvarkit;isHg19;isHg38;isBlank;hasFeature;moduleLoad;getGnomadGenomePath;getGnomadExomePath} from '../../modules/utils/functions.nf'
include {GO_DOWNLOAD_OBO_01} from '../../modules/go/go.download.obo.01.nf'
include {GOA_DOWNLOAD_01} from '../../modules/goa/goa.download.01.nf'
include {CONCAT_FILES_01} from '../../modules/utils/concat.files.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {DOWNLOAD_GFF3_01} from '../../modules/gff3/download.gff3.01.nf'



workflow GO_TO_BED_01 {
	take:
		meta
		reference
		gff3
		whiteListedGeneList
		blackListedGeneList
	main:
		version_ch = Channel.empty()

		obo_ch = GO_DOWNLOAD_OBO_01([:])
		version_ch = version_ch.mix(obo_ch.version)

		goa_ch = GOA_DOWNLOAD_01([:], reference)
		version_ch = version_ch.mix(goa_ch.version)

		if(gff3.name.equals("NO_FILE")) {
			gff3_ch = DOWNLOAD_GFF3_01(meta.plus(with_tabix:"true"),reference)
			version_ch = version_ch.mix(gff3_ch.version)
			gff = gff3_ch.gff3
			}
		else
			{
			gff = gff3
			}

		bed4genes_ch = BED_FOR_GENES(meta, reference, gff, obo_ch.obo, goa_ch.gaf, whiteListedGeneList, blackListedGeneList)
                version_ch = version_ch.mix(bed4genes_ch.version)

		version_ch = MERGE_VERSION(meta, "go2bed", "GeneOntology To Bed", version_ch.collect())
	emit:
		version = version_ch
		obo = obo_ch.obo
		goa = goa_ch.gaf
		bed = bed4genes_ch.bed
		index = bed4genes_ch.bed_index
		terms = bed4genes_ch.terms
	}


/* gene symbols as bed file sorted on gene symbol */
process BED_FOR_GENES {
afterScript "rm -rf TMP"
memory "10g"
input:
	val(meta)
	val(reference)
	path(gff3)
      	path(obo)
	path(goa)
	val(whiteListedGeneList)
	val(blackListedGeneList)
output:
       	path("go.bed.gz"),emit:bed
       	path("go.bed.gz.tbi"),emit:bed_index
       	path("go.flag.bed.gz"),emit:flag_bed
	path("go.tsv"),emit:output //output for bcftools annotate
	path("go.flag.tsv"),emit:flag_output //output for bcftools annotate, flag only
	path("go.terms.tsv"),emit:terms
	path("version.xml"),emit:version
script:
	def whatis="Genes for Gene Ontology terms ${meta.goTerms} . Exclude terms = ${meta.excludeGoTerms}"
"""
hostname 1>&2
${moduleLoad("jvarkit htslib datamash bedtools")}
set -o pipefail
mkdir -p TMP

export LC_ALL=C

echo '${meta.goTerms}' | tr " ,|;" "\\n" | sort |uniq | grep -v '^\$' > TMP/include.terms

if ${!meta.excludeGoTerms.trim().isEmpty()} ; then
	echo '${meta.excludeGoTerms}' | tr " ,|;" "\\n" | sort |uniq | grep -v '^\$' > TMP/exclude.terms
fi

# extract genes, sort by gene name
java  -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/gtf2bed.jar --columns "gtf.feature,Name,gene_id" -R "${reference}" "${gff3}" |\
awk -F '\t' '(\$4=="gene" && \$5!=".")' |\
cut -f1,2,3,5,6 |\
LC_ALL=C sort -T TMP -t '\t' -k4,4 | uniq > TMP/genes.symbols.bed

# extract go terms ACN/NAME/DESCRIPTION/DIVISION
java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/goutils.jar \
	--action dump_table \
	--accession-file "TMP/include.terms" \
	-go "${obo}"  |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 > TMP/go.terms.tsv


if test -s TMP/exclude.terms ; then

	# extract go terms ACN/NAME/DESCRIPTION/DIVISION
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/goutils.jar \
		--action dump_table \
		--accession-file "TMP/exclude.terms" \
		-go "${obo}"  |
		cut -f1 | LC_ALL=C sort -T TMP -t '\t' -k1,1 | uniq  > TMP/x.terms.tsv

	LC_ALL=C join -t '\t' -1 1 -2 1 -o '1.1,1.2,1.3,1.4' -v1 TMP/go.terms.tsv TMP/x.terms.tsv > TMP/jeter1.txt

	mv TMP/jeter1.txt TMP/go.terms.tsv

fi

# sort GOA on gene name
grep -v '^!'  '${goa}' | awk -F '\t' '!(\$4 ~ /^NOT\\|/)' |\
	cut -f3,5 | LC_ALL=C sort -T TMP -t '\t' -k2,2 | uniq > TMP/sn_acn.tsv

# join terms an GOA
LC_ALL=C join -t '\t' -1 1 -2 2  -o '1.1,1.2,2.1' TMP/go.terms.tsv  TMP/sn_acn.tsv |\
	awk -F '\t' '{OFS="\t";gsub(/[^A-Za-z0-9]+/,"_",\$2);print}' |\
	LC_ALL=C sort -T TMP -t '\t' -k3,3  > TMP/gene_go.tsv

# join on BED
LC_ALL=C join -t '\t' -1 4 -2 3 -o '1.1,1.2,1.3,1.4,1.5,2.1,2.2' TMP/genes.symbols.bed  TMP/gene_go.tsv |\
	sed 's/\t/|/6' |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 |\
	LC_ALL=C datamash groupby 1,2,3,4,5  unique 6 >  TMP/jeter.bed

# white list
if ${!file(whiteListedGeneList).name.equals("NO_FILE")} ; then

	cat "${whiteListedGeneList}" | sort -T TMP | uniq > TMP/jeter.a
	# gene already in bed
	cut -f 4 TMP/jeter.bed |  sort -T TMP | uniq > TMP/jeter.b
	# unique to white list
	comm -23 TMP/jeter.a TMP/jeter.b > TMP/jeter.c
	
	if test -s TMP/jeter.c ; then
		
		join -t '\t' -1 4 -2 1 TMP/genes.symbols.bed TMP/jeter.c -o '1.1,1.2,1.3,1.4,1.5' |\
			sed 's/\$/\t./' >> TMP/jeter.bed

		LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n TMP/jeter.bed > TMP/jeter2.bed
		mv TMP/jeter2.bed TMP/jeter.bed
	fi
fi


# blacklisted list
if ${!file(blackListedGeneList).name.equals("NO_FILE")} ; then

	cat "${blackListedGeneList}" | sort -T TMP | uniq > TMP/jeter.a
		
	sort -t '\t' -T TMP -k4,4 TMP/jeter.bed |\
		join -t '\t' -1 4 -2 1 -v 1  -o '1.1,1.2,1.3,1.4,1.5,1.6'  - TMP/jeter.a  > TMP/jeter2.bed
	mv TMP/jeter2.bed TMP/jeter.bed

fi

LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n TMP/jeter.bed > TMP/jeter2.bed
mv TMP/jeter2.bed TMP/jeter.bed

bgzip TMP/jeter.bed

mv TMP/jeter.bed.gz go.bed.gz

tabix -p bed go.bed.gz

mv TMP/go.terms.tsv ./

echo '##INFO=<ID=GO,Number=.,Type=String,Description="${whatis}">' > go.header

echo "\${PWD}/go.bed.gz\tCHROM,FROM,TO,-,-,GO\t\${PWD}/go.header" > go.tsv


gunzip -c go.bed.gz | cut -f1,2,3 | bedtools merge | sed 's/\$/\t1/' | bgzip > go.flag.bed.gz
tabix -p bed go.flag.bed.gz

echo '##INFO=<ID=GO,Number=1,Type=Flag,Description="${whatis}">' > go.flag.header
echo "\${PWD}/go.flag.bed.gz\tCHROM,FROM,TO,GO\t\${PWD}/go.flag.header" > go.flag.tsv


#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">${whatis}</entry>
        <entry key="obo">${obo}</entry>
        <entry key="goa">${goa}</entry>
        <entry key="gff3">${gff3}</entry>
        <entry key="include.go">${meta.goTerms}</entry>
        <entry key="exclude.go">${meta.excludeGoTerms}</entry>
        <entry key="versions">${getVersionCmd("jvarkit/goutils jvarkit/gtf2bed tabix datamash awk")}</entry>
</properties>
EOF
"""
}
