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

include {moduleLoad;assertFileExists;isBlank} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {DOWNLOAD_SOMALIER} from '../../modules/somalier/somalier.download.nf'
include {SOMALIER_DOWNLOAD_SITES} from '../../modules/somalier/somalier.download.sites.nf'
include {SAMTOOLS_SAMPLES} from '../samtools/samtools.samples.03.nf'




workflow SOMALIER_BAMS_01 {
	take:
		genome_ch // fasta,fai,dict
		samplesheet // sample,bam,bai
		pedigree // pedigree for somalier
		user_sites //file or no file
	main:
		version_ch = Channel.empty()

		exe_ch = DOWNLOAD_SOMALIER()
		version_ch = version_ch.mix(exe_ch.version)


		if(user_sites.name.equals("NO_FILE")) {
			sites_ch = SOMALIER_DOWNLOAD_SITES(genome_ch)
			version_ch = version_ch.mix(sites_ch.version)
			sites_vcf= sites_ch.output
			}
		else {
			sites_vcf = Channel.of(user_sites).
				map{[file(it),file(it+".tbi")]}
			}

		ch3 = EXTRACT_BAM(genome_ch, exe_ch.output , samplesheet.combine(sites_vcf))
		version_ch = version_ch.mix(ch3.version)
	
		somalier_ch = RELATE_SOMALIER(genome_ch, exe_ch.output, ch3.output.collect(), pedigree)
		version_ch = version_ch.mix(somalier_ch.version)

			
		//version_ch = MERGE_VERSION( "somalier",version_ch.collect())

	emit:
		output = somalier_ch.output
		version = version_ch
		zip = somalier_ch.zip
		qc = somalier_ch.qc

	}


process EXTRACT_BAM {
	tag "${sample}"
	label "process_quick"
	memory '2g'
	input:
		tuple path(fasta),path(fai),path(dict)
		path(somalier)
		tuple val(sample),path(bam),path(bai),path(sites),path(sites_idx)
	output:
		path("extracted/${sample}.somalier"),emit:output
		path("version.xml"),emit:version
	script:
	"""
	hostname 1>&2
	mkdir -p extracted
	./${somalier} extract -d extracted --sites "${sites}" -f "${fasta}" "${bam}"
	
	test -s "extracted/${sample}.somalier"

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">somalier extract</entry>
        <entry key="sample">${sample}</entry>
        <entry key="bam">${bam}</entry>
</properties>
EOF
	"""
	}


process RELATE_SOMALIER {
afterScript "rm -rf extracted TMP"
input:
	tuple path(fasta),path(fai),path(dict)
	path(somalier)
	path("EXTRACT/*")
	path(pedigree)
output:
	path("${params.prefix?:""}somalier.bams/*"),emit:output
	path("${params.prefix?:""}somalier.bams.zip"),emit:zip
	path("${params.prefix?:""}somalier_mqc.html"),emit:qc
	path("version.xml"),emit:version
script:
	def max_rows_html = 50
	def prefix=params.prefix?:""
"""
hostname 1>&2
set -x

mkdir -p TMP
mkdir -p "${prefix}somalier.bams"

./${somalier} relate --output-prefix=${prefix}somalier.bams/${prefix}bams \\
	${pedigree.name.equals("NO_FILE")?"":"-p '${pedigree}'"} \\
	EXTRACT/*.somalier

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
        <entry key="pedigree">${pedigree}</entry>
</properties>
EOF
"""
}
