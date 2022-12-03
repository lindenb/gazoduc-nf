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
include {GO_DOWNLOAD_OBO_01} from '../../modules/go/go.download.obo.01.nf'
include {GOA_DOWNLOAD_01} from '../../modules/goa/goa.download.01.nf'
include {CONCAT_FILES_01} from '../../modules/utils/concat.files.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'


workflow GO_TO_BED_01 {
	take:
		meta
		reference
		gtf
	main:
		version_ch = Channel.empty()

		obo_ch = GO_DOWNLOAD_OBO_01(meta)
		version_ch = version_ch.mix(obo_ch.version)

		goa_ch = GOA_DOWNLOAD_01(meta, reference)
		version_ch = version_ch.mix(goa_ch.version)

		if(gtf.name.equals("NO_FILE")) {
			sorted_bed = bed
			}
		else
			{
			gtf1 = gtf
			}

		bed4genes_ch = BED_FOR_GENES(meta, reference, gtf1)
                version_ch = version_ch.mix(bed4genes_ch.version)

		version_ch = MERGE_VERSION(meta, "go2bed", "GeneOntology To Bed", version_ch.collect())
	emit:
		version = version_ch
		obo = obo_ch.obo
		goa = goa_ch.gaf
		bed = 
	}


/* gene symbols as bed file sorted on gene symbol */
process BED_FOR_GENES {
afterScript "rm -rf TMP"
memory "10g"
input:
	val(meta)
	val(reference)
	path(gtf)
      	path(obo)
	path(goa)
output:
       	path("genes.symbols.bed"),emit:bed
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bedtools jvarkit htslib")}
set -o pipefail
mkdir -p TMP


echo "${meta.goTerms}" | tr " ,|;" "\\n" | sort |uniq | grep -v '^\$' > TMP/jeter.terms
test -s TMP/jeter.terms


java -jar \${JVARKIT_DIST}/gtf2bed.jar --columns "gtf.feature,gene_name,gene_id" -R "${reference}" "${gtf}" |\
awk -F '\t' '(\$4=="gene")' |\
cut -f1,2,3,5,6 |\
sort -T TMP -t '\t' -k4,4| uniq > TMP/genes.symbols.bed


java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/goutils.jar \
	--action dump_table \
	--accession-file "jeter.terms" \
	-go "${obo}"  | grep -F -f jeter.terms > go.terms.tsv

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
