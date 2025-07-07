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
include {k1_signature} from '../../../modules/utils/k1.nf'



workflow VISTA {
	take:
		meta
		fasta
		fai
		dict
		vcfs /* meta, vcf,vcf_index */
	main:
		versions = Channel.empty()
		
		DOWNLOAD(fasta,fai,dict)
		versions = versions.mix(DOWNLOAD.out.versions)

		ANNOTATE(DOWNLOAD.out.bed, DOWNLOAD.out.tbi,DOWNLOAD.out.header,vcfs)
		versions = versions.mix(ANNOTATE.out.versions)
	emit:
		vcf = ANNOTATE.out.vcf
		versions
		doc = DOWNLOAD.out.doc
}

process DOWNLOAD{
tag "${meta1.id?:fasta.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
output:
	tuple val(meta1),path("*.bed.gz"),emit:bed
	tuple val(meta1),path("*.bed.gz.tbi"),emit:tbi
	tuple val(meta1),path("*.header"),emit:header
	tuple val(meta1),path("*.md"),emit:doc
	tuple val(meta1),path("versions.yml"),emit:versions
script:
    def k1 = k1_signature()
   	def TAG = "VISTA"
	def whatis="VISTA enhancers"
	def base = "https://hgdownload.cse.ucsc.edu/gbdb"
"""
set -o pipefail
hostname 1>&2
mkdir -p TMP/CACHE

cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg38}\t${base}/hg38/vistaEnhancers/vistaEnhancers.bb
1:${k1.hg19}\t${base}/hg19/vistaEnhancers/vistaEnhancers.bb
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sed 's/^chr//' | sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv
join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv | sort | uniq > TMP/jeter.url

test -s TMP/jeter.url

wget -O TMP/jeter.bb `cat TMP/jeter.url`

bigBedToBed -udcDir=TMP/CACHE TMP/jeter.bb stdout |\\
        cut -f1,2,3,4 |\\
        jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP bedrenamechr -R ${fasta} --column 1 --convert SKIP  |\\
        LC_ALL=C sort -t '\t' -k1,1 -k2,2n -T TMP |\\
        bgzip > TMP/${TAG}.bed.gz

tabix -p bed -f TMP/${TAG}.bed.gz

mv TMP/${TAG}.bed.gz ./
mv TMP/${TAG}.bed.gz.tbi ./

echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="${whatis}">' > ${TAG}.header

cat << EOF > ${TAG}.md
Vista enhancer.
EOF

cat << END_VERSIONS > versions.yml
"${task.process}":
	TODO: "TODO"
END_VERSIONS
"""
}

process ANNOTATE {
tag "${meta.id?:vcf.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(tabix)
	tuple val(meta2),path(tbi)
	tuple val(meta3),path(header)
	tuple val(meta),path(vcf),path(vcf_idx)
output:
    tuple val(meta),path("*.bcf"),path("*.csi"),emit:vcf
	path("versions.yml"),emit:versions
script:
    def TAG = "VISTA"
"""
hostname 1>&2
mkdir -p TMP OUTPUT

bcftools annotate \\
    --threads ${task.cpus} \\
    -a "${tabix}" \\
    -h "${header}" \\
    -c "CHROM,FROM,TO,${TAG}" \\
    --merge-logic '${TAG}:unique' -O b -o TMP/${TAG}.${vcf.getBaseName()}.bcf '${vcf}'
    
bcftools index \\
    --threads ${task.cpus} \\
    --force TMP/${TAG}.${vcf.getBaseName()}.bcf

mv TMP/*.bcf ./
mv TMP/*.bcf.csi ./


cat << END_VERSIONS > versions.yml
"${task.process}":
	TODO: "TODO"
END_VERSIONS
"""
}
