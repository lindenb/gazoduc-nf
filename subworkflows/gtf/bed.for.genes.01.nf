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
include { getKeyValue; moduleLoad; getBoolean; parseBoolean} from '../../modules/utils/functions.nf'
include {COMPILE_GTF_TO_BED} from '../../modules/gtf/gtf2bed.compile.01.nf' 
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {GTF_TO_GENES_BED_01} from '../../modules/gtf/gtf2bed.genes.01.gtf'

workflow BED_FOR_GENES_01 {
	take:
		meta
		reference
		user_genes
	main:
		version_ch = Channel.empty()
		
		jar_ch = COMPILE_GTF_TO_BED([:])
		version_ch= version_ch.mix(jar_ch.version)

		gene_bed_ch = GTF_TO_GENES_BED_01(meta, reference, jar_ch.jar )
		version_ch= version_ch.mix(gene_bed_ch.version)

		join_ch = JOIN(meta, gene_bed_ch.bed, user_genes )
		version_ch= version_ch.mix(join_ch.version)

		version_merge = MERGE_VERSION(meta,"bed4gene","Bed from user genes",version_ch.collect())
		
	emit:
		version = version_merge.version
		bed = join_ch.bed
	}


process JOIN {
executor "local"
afterScript "rm -rf TMP"
input:
	val(meta)
	path(genesbed)
	path(usergenes)
output:
	path("user_genes.bed"),emit:bed
	path("version.xml"),emit:version
script:
	def check_unpaired = parseBoolean(getKeyValue(meta,"check_unpaired","true"))
	def with_header = parseBoolean(getKeyValue(meta,"with_header","false"))
"""
	hostname 1>&2
	set -o pipefail
	mkdir TMP
	test -s "${usergenes}"

	grep -v '^#' "${genesbed}" | LC_ALL=C sort -T . -t '\t' -k6,6 > TMP/a
	cut -f 1 "${usergenes}" | grep -v '^\$' | LC_ALL=C sort -T . -t '\t' -k1,1 > TMP/b

	LC_ALL=C join -t '\t' -o '1.1,1.2,1.3,1.4,1.5,1.6' -1 6 -2 1 TMP/a TMP/b |\
		LC_ALL=C sort -t '\t' -T TMP -k1,1 -k2,2n > TMP/jeter.bed

	cut -f 6 TMP/jeter.bed | sort | uniq -d > TMP/dups.txt
	echo "Duplicate names:" 1>&2
	cat TMP/dups.txt 1>&2
	grep -F -f TMP/dups.txt TMP/jeter.bed 1>&2 || true
	test ! -s  TMP/dups.txt

	if [ ! -z "${check_unpaired?"Y":""}" ] ; then	
		LC_ALL=C join -t '\t' -v 2  -1 6 -2 1 TMP/a TMP/b > TMP/unpaired.txt
		echo "#UNPAIRED USER GENES:" 1>&2
		cat TMP/unpaired.txt 1>&1
		test ! -s TMP/unpaired.txt
	fi

	if [ ! -z "${with_header?"Y":""}" ] ; then
		echo -e "#chrom\tstart\tend\tgene_type\tgene_id\tgene_name" > TMP/jeter2.bed
		cat TMP/jeter.bed >> TMP/jeter2.bed
		mv TMP/jeter2.bed TMP/jeter.bed
	fi
	
	mv TMP/jeter.bed user_genes.bed

	######################################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">join list of gene and user data</entry>
		<entry key="genes">${genesbed}</entry>
		<entry key="user.genes">${usergenes}</entry>
		<entry key="user.genes.count">\$(wc -l < ${usergenes})</entry>
	</properties>
	EOF
"""
}
