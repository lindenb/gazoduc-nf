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

include {moduleLoad;assertFileExists;isBlank} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {DOWNLOAD_SOMALIER} from '../../modules/somalier/somalier.download.nf'
include {SOMALIER_DOWNLOAD_SITES} from '../../modules/somalier/somalier.download.sites.nf'
include {SAMTOOLS_SAMPLES} from '../samtools/samtools.samples.03.nf'




workflow SOMALIER_BAMS_01 {
	take:
		meta
		genomeId
		bams
		pedigree
		user_sites
	main:
		version_ch = Channel.empty()

		ch1 = SAMTOOLS_SAMPLES([:], bams)
		version_ch = version_ch.mix(ch1.version)

		ch2 = ch1.rows.filter{T->T.genomeId.equals(genomeId)}

		som_ch = SOMALIER_BAMS_02([:], genomeId, ch2, pedigree, user_sites)
		version_ch = version_ch.mix(som_ch.version)

		version_ch = MERGE_VERSION( "somalier",version_ch.collect())
	emit:
		version = version_ch
		zip = som_ch.zip
		output = som_ch.output
		qc = som_ch.qc /** input for multiqc */
	}

workflow SOMALIER_BAMS_02 {
	take:
		meta
		genomeId
		rows
		pedigree
		user_sites
	main:
		version_ch = Channel.empty()

		exe_ch = DOWNLOAD_SOMALIER([:])
		version_ch = version_ch.mix(exe_ch.version)

		if(user_sites.name.equals("NO_FILE")) {
			sites_ch = SOMALIER_DOWNLOAD_SITES([:],genomeId)
			version_ch = version_ch.mix(sites_ch.version)
			sites_vcf= sites_ch.vcf
			}
		else {
			sites_vcf = user_sites
			}

		ch3 = EXTRACT_BAM([:], genomeId, exe_ch.executable , sites_vcf , rows)
		version_ch = version_ch.mix(ch3.version)
			
		somalier_ch = RELATE_SOMALIER([:], genomeId, exe_ch.executable,ch3.output.collect(), pedigree)
		version_ch = version_ch.mix(somalier_ch.version)
	
		version_ch = MERGE_VERSION( "somalier",version_ch.collect())
	emit:
		output = somalier_ch.output
		version = version_ch
		zip = somalier_ch.zip
		qc = somalier_ch.qc
	}


process EXTRACT_BAM {
	tag "${row.sample}"
	memory '2g'
	input:
		val(meta)
		val(genomeId)
		val(somalier)
		val(sites)
		val(row)
	output:
		path("extracted/${row.sample}.somalier"),emit:output
		path("version.xml"),emit:version
	script:
		def genome = params.genomes[genomeId]
	"""
	hostname 1>&2
	mkdir -p extracted
	${somalier} extract -d extracted --sites "${sites}" -f "${genome.fasta}" "${row.bam}"
	
	test -s "extracted/${row.sample}.somalier"

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">somalier extract</entry>
        <entry key="sample">${row.sample}</entry>
        <entry key="bam">${row.bam}</entry>
</properties>
EOF
	"""
	}


process RELATE_SOMALIER {
tag "N=${L.size()}"
afterScript "rm -rf extracted TMP"
input:
	val(meta)
	val(genomeId)
	val(somalier)
	val(L)
	path(pedigree)
output:
	path("${params.prefix?:""}somalier.bams/*"),emit:output
	path("${params.prefix?:""}somalier.bams.zip"),emit:zip
	path("${params.prefix?:""}somalier_mqc.html"),emit:qc
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def prefix = params.prefix?:""
	def max_rows_html = 50
"""
hostname 1>&2
set -x

mkdir -p TMP
mkdir -p "${prefix}somalier.bams"

${somalier} relate --output-prefix=${prefix}somalier.bams/${prefix}bams \
	${pedigree.name.equals("NO_FILE")?"":"-p '${pedigree}'"} \
	${L.join(" ")}

# may not exist
touch "${prefix}somalier.bams/${prefix}bams.groups.tsv"
zip -9 -r "${prefix}somalier.bams.zip" "${prefix}somalier.bams"


cat << EOF > TMP/jeter.html
<!--
id: '${prefix}'
section_name: 'Somalier'
description: 'Somalier: relatedness among samples from extracted, genotype-like information'
-->
<div>
<table class="table">
<caption>${max_rows_html} higher relatedness</caption>
<thead>
	<th>sample1</th>
	<th>sample2</th>
	<th>relatedness</th>
</thead>
<tbody>
EOF

cut -f1,2,3  "${prefix}somalier.bams"/*.pairs.tsv |\
	tail -n+2 |\
	LC_ALL=C sort -T TMP -t '\t' -k3,3gr |\
	head -n '${max_rows_html}' |\
	awk -F '\t' '{printf("<tr><td>%s</td><td>%s</td><td>%s</td></tr>\\n",\$1,\$2,\$3);}' >> TMP/jeter.html
cat << EOF >> TMP/jeter.html
</tbody>
</table>
<br/>

<table class="table">
<caption>mean relatedness by sample</caption>
<thead>
        <th>sample</th>
        <th>AVG(relatedness)</th>
</thead>
<tbody>
EOF

awk '(NR>1) {P[\$1]+=1.0*(\$3);P[\$2]+=1.0*(\$3);C[\$1]++;C[\$2]++;} END{for(S in P) printf("%s\t%f\\n",S,P[S]/C[S]);}' "${prefix}somalier.bams"/*.pairs.tsv  |\
	LC_ALL=C sort -T TMP -t \$'\\t' -k2,2gr |\
	awk -F '\t' '{printf("<tr><td>%s</td><td>%s</td></tr>\\n",\$1,\$2);}'  >> TMP/jeter.html

cat << EOF >> TMP/jeter.html
</tbody>
</table>
</div>
EOF

mv -v  TMP/jeter.html "${prefix}somalier_mqc.html"

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">run somalier on bams</entry>
        <entry key="bams.count">${L.size()}</entry>
        <entry key="pedigree">${pedigree}</entry>
</properties>
EOF
"""
}
