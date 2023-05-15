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
include {getVersionCmd;moduleLoad;parseBoolean} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {GTF_TO_GENES_BED_02} from '../../modules/gtf/gtf2bed.genes.02.gtf'


/** build a BED for bcftools annotate with jensenlab diseases */
workflow JENSENLAB_DISEASES_01 {
	take:
		meta
		reference
		gtf
	main:

		version_ch  = Channel.empty()

		ch1 = GTF_TO_GENES_BED_02(["with_header":true], reference, gtf)
		version_ch = version_ch.mix(ch1.version)
	
		ch2= BUILD([:], reference, ch1.bed)
		version_ch = version_ch.mix(ch2.version)

		version_ch = MERGE_VERSION(meta, "jensenlab.diseases", "jensenlab.diseases", version_ch.collect())

	emit:
		output = ch2.output
		version = version_ch
		bed = ch2.bed
		header = ch2.header
	}


process BUILD {
	tag "${bed}"
	input:
		val(meta)
		val(reference)
		val(bed)
	output:
		path("diseases.tsv"),emit:output
		path("diseases.bed.gz"),emit:bed
		path("diseases.bed.gz.tbi"),emit:index
		path("diseases.header"),emit:header
		path("version.xml"),emit:version
	script:
		def url = "https://download.jensenlab.org/human_disease_textmining_filtered.tsv"
	"""
	hostname 1>&2
	${moduleLoad("htslib")}
	set -o pipefail
	

	mkdir -p TMP
	gunzip -c "${bed}" |\
		awk -F '\t' 'BEGIN{C=-1;} /^#/ {for(i=3;i<=NF && C<0;i++){if(\$i=="gene_name") {C=i;break;}}} {if(C>0) printf("%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,\$C);}' |\
		LC_ALL=C sort -T TMP -t '\t' -k4,4 > TMP/jeter.a

	test -s TMP/jeter.a

	 wget  -O - "${url}" |\
		awk -F '\t' 'BEGIN {printf("gene_name\tdoid\tdisease\\n");} {D=\$3;gsub(/[\\:]/,"_",D);M=\$4;gsub(/[^A-Za-z0-9]+/,"_",M);printf("%s\t%s\t%s\\n",\$2,D,M);}' |\
		LC_ALL=C sort -T TMP -t '\t' -k1,1 > TMP/jeter.b

	test -s TMP/jeter.b


	LC_ALL=C join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,1.4,2.2,2.3' TMP/jeter.a TMP/jeter.b |\
		LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n | uniq > TMP/jeter.bed

	bgzip TMP/jeter.bed
        tabix --comment '#' -f -p bed TMP/jeter.bed.gz

	mv TMP/jeter.bed.gz diseases.bed.gz
	mv TMP/jeter.bed.gz.tbi diseases.bed.gz.tbi

	echo '##INFO=<ID=DISEASES,Number=.,Type=String,Description="DISEASES is a weekly updated web resource that integrates evidence on disease-gene associations from automatic text mining, manually curated literature, cancer mutation data, and genome-wide association studies. ${url}">' > diseases.header
	
	echo "\${PWD}/diseases.bed.gz\tCHROM,FROM,TO,GENE_NAME,DOID,DISEASES\t\${PWD}/diseases.header" > diseases.tsv

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">join diseaes and bed gene</entry>
        <entry key="url"><a>${url}</a></entry>
</properties>
EOF
	"""
	}

