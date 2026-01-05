/*

Copyright (c) 2026 Pierre Lindenbaum

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


workflow REMAP {
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

		ANNOTATE(DOWNLOAD.out.output,vcfs)
        versions = versions.mix(ANNOTATE.out.versions)
	emit:
		vcf = ANNOTATE.out.vcf
		versions
	}


process DOWNLOAD {
tag "${meta.id?:fasta.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
output:
	tuple val(meta1),
        path("*.bed.gz"), path("*.bed.gz.tbi"), path("*.header"),emit:output
script:
    def k1 = k1_signature()
    def TAG = "REMAP"
    def WHATIZ = "ReMap is a large scale integrative analysis of DNA-binding experiments for Homo sapiens, Mus musculus, Drosophila melanogaster and Arabidopsis thaliana transcriptional regulators. The catalogues are the results of the manual curation of ChIP-seq, ChIP-exo, DAP-seq from public sources (GEO, ENCODE, ENA)."
    def base = "https://remap.univ-amu.fr/storage/remap2022/"
"""
hostname 1>&2
mkdir -p TMP
set -o pipefail
mkdir -p TMP/CACHE


cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg38}\t${base}/hg38/MACS2/remap2022_crm_macs2_hg38_v1_0.bed.gz
1:${k1.hg19}\t${base}/hg19/MACS2/remap2022_crm_macs2_hg19_v1_0.bed.gz
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' |\\
    sed 's/^chr//' |\\
    sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv
join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv |\\
sort | uniq > TMP/jeter.url

test -s TMP/jeter.url

wget -O - `cat TMP/jeter.url`  |\\
        gunzip -c |\\
        cut -f1,2,3 |\\
        jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP bedrenamechr -R "${fasta}" --column 1 --convert SKIP  |\\
        LC_ALL=C sort --buffer-size=${task.memory.mega}M -t '\t' -k1,1 -k2,2n -T TMP |\\
        bedtools merge |\\
	sed 's/\$/\t1/' |\\
        bgzip > TMP/${TAG}.bed.gz


tabix -p bed -f TMP/${TAG}.bed.gz

mv TMP/${TAG}.bed.gz ./
mv TMP/${TAG}.bed.gz.tbi ./


echo '##INFO=<ID=${TAG},Number=0,Type=Flag,Description="${WHATIZ}">' > ${TAG}.header

cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "\${URL}"
END_VERSIONS
"""
}

process ANNOTATE {
tag "${meta.id?:vcf.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(tabix),path(tbi),path(header)
	tuple val(meta),path(vcf),path(vcf_idx)
output:
    tuple val(meta),path("*.bcf"),path("*.csi"),emit:output
script:
    def TAG = "REMAP"
    def prefix = task.ext.prefix?:vcf.baseName+".remap"
"""
hostname 1>&2
mkdir -p TMP OUTPUT

bcftools annotate \\
    --write-index \\
    --threads ${task.cpus} \\
    -a "${tabix}" \\
    -h "${header}" \\
    -c "CHROM,FROM,TO,${TAG}" \\
    -O b \\
    -o TMP/${prefix}.bcf '${vcf}'

mv TMP/*.bcf ./
mv TMP/*.bcf.csi ./

cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
}
