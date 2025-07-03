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

//include {DOWNLOAD_SOMALIER} from '../../modules/somalier/somalier.download.nf'
include {SOMALIER_DOWNLOAD_SITES} from '../../../modules/somalier/download.sites/main.nf'
//include {SAMTOOLS_SAMPLES} from '../samtools/samtools.samples.03.nf'




workflow SOMALIER_BAMS {
	take:
		meta
		fasta
		fai
		dict
		samplesheet // sample,bam,bai
		pedigree // pedigree for somalier
		user_sites //file or no file
	main:
		version_ch = Channel.empty()



		if(user_sites[1]) {
			sites_vcf = Channel.of(user_sites).
				map{[file(it),file(it+".tbi")]}.collect()
			}
		else
			{
			SOMALIER_DOWNLOAD_SITES(fasta,fai,dict)
			version_ch = version_ch.mix(SOMALIER_DOWNLOAD_SITES.out.versions)
			sites_vcf= SOMALIER_DOWNLOAD_SITES.out.vcf
			}

		EXTRACT_BAM(fasta,fai,dict, sites_vcf, samplesheet)
		version_ch = version_ch.mix(EXTRACT_BAM.out.versions.first())
	
		RELATE_SOMALIER(
			fasta,fai,dict,
			EXTRACT_BAM.out.output.collect(),
			pedigree
			)
		version_ch = version_ch.mix(RELATE_SOMALIER.out.versions)

	emit:
		output = RELATE_SOMALIER.out.output
		versions = version_ch
		zip = somalier_ch.zip
		qc = somalier_ch.qc

	}


process EXTRACT_BAM {
	tag "${meta.id?:bam.name}"
	label "process_quick"
	conda "${moduleDir}/../../conda/somalier.yml"
	memory '2g'
	input:
		tuple val(meta1),path(fasta)
		tuple val(meta2),path(fai)
		tuple val(meta3),path(dict)
		tuple val(meta4),path(sites),path(sites_idx)
		tuple val(meta),path(bam),path(bai)
	output:
		path("extracted/*.somalier"),emit:output
		path("versions.yml"),emit:versions
	script:
	"""
	hostname 1>&2
	mkdir -p extracted
	somalier extract -d extracted --sites "${sites}" -f "${fasta}" "${bam}"
	
	test -s "extracted/${meta.id}.somalier"


cat << EOF > versions.yml
${task.process}:
    somalier: todo
EOF
	"""
	}


process RELATE_SOMALIER {
afterScript "rm -rf extracted TMP"
conda  "${moduleDir}/../../conda/somalier.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	path("EXTRACT/*")
	tuple val(meta4),path(pedigree)
output:
	path("somalier.bams/*"),emit:output
	path("somalier.bams.zip"),emit:zip
	path("somalier_mqc.html"),emit:qc
	path("versions.yml"),emit:versions
script:
	def max_rows_html = 50
	def prefix=task.ext.prefix?:""
"""
hostname 1>&2
set -x

mkdir -p TMP
mkdir -p "${prefix}somalier.bams"

somalier relate --output-prefix=${prefix}somalier.bams/${prefix}bams \\
	${pedigree?"-p '${pedigree}'":""} \\
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

cat << EOF > versions.yml
${task.process}:
    somalier: todo
EOF
"""
}
