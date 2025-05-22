/*

Copyright (c) 2025 Pierre Lindenbaum

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


workflow ANNOTATE_REMAP {
	take:
		fasta
		fai
		dict
		vcfs /* meta, vcf,vcf_index */
	main:
		source_ch = DOWNLOAD(fasta,fai,dict)
		annotate_ch = ANNOTATE(source_ch.bed, source_ch.tbi,source_ch.header,vcfs)
	emit:
		output = annotate_ch.output
		doc = source_ch.doc
	}


process DOWNLOAD {
tag "${fasta.name}"
afterScript "rm -rf TMP"
label "process_quick"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
        path(fasta)
        path(fai)
        path(dict)
output:
	path("*.bed.gz"),emit:bed
	path("*.bed.gz.tbi"),emit:tbi
	path("*.header"),emit:header
	path("*.md"),emit:doc
script:
    def k1 = k1_signature()
    def TAG = "REMAP"
    def WHATIZ = "ReMap is a large scale integrative analysis of DNA-binding experiments for Homo sapiens, Mus musculus, Drosophila melanogaster and Arabidopsis thaliana transcriptional regulators. The catalogues are the results of the manual curation of ChIP-seq, ChIP-exo, DAP-seq from public sources (GEO, ENCODE, ENA)."
"""
hostname 1>&2
mkdir -p TMP
set -o pipefail
mkdir -p TMP/CACHE


cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg38}\thttps://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/remap2022_crm_macs2_hg38_v1_0.bed.gz
1:${k1.hg19}\thttps://remap.univ-amu.fr/storage/remap2022/hg19/MACS2/remap2022_crm_macs2_hg19_v1_0.bed.gz
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sed 's/^chr//' | sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv
join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv | sort | uniq > TMP/jeter.url

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

cat << __EOF__ > ${TAG}.md
${WHATIZ}
__EOF__

"""
}

process ANNOTATE {
tag "${vcf.name}"
afterScript "rm -rf TMP"
label "process_quick"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	path(tabix)
	path(tbi)
	path(header)
	tuple val(meta),path(vcf),path(vcf_idx)
output:
    tuple val(meta),path("*.bcf"),path("*.csi"),emit:output
script:
    def TAG = "REMAP"
"""
hostname 1>&2
mkdir -p TMP OUTPUT

bcftools annotate \\
    --threads ${task.cpus} \\
    -a "${tabix}" \\
    -h "${header}" \\
    -c "CHROM,FROM,TO,${TAG}" \\
    -O b -o TMP/${TAG}.${vcf.getSimpleName()}.bcf '${vcf}'
    
bcftools index \\
    --threads ${task.cpus} \\
    --force TMP/${TAG}.${vcf.getSimpleName()}.bcf

mv TMP/*.bcf ./
mv TMP/*.bcf.csi ./
"""
}
