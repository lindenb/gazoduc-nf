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



workflow ANNOTATE_HEART_ENHANCER_DICKEL {
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

process DOWNLOAD{
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
   	def TAG = "ENHANCER_DICKEL"
	def url = "https://portal.nersc.gov/dna/RD/heart/literature_metaAnalysis_trackHub/hg19/dickel_et_al_scores/output.scoreAll.bb"
	def whatis="Heart enhancer predictions (Dickel et al.) ${url}"
"""
hostname 1>&2
mkdir -p TMP/CACHE
set -o pipefail

wget -O TMP/jeter.bb "${url}"

cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg38}\thg38
1:${k1.hg19}\thg39
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sed 's/^chr//' | sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv
join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv | sort | uniq > TMP/jeter.build

test -s TMP/jeter.build

if grep -F -w hg38 TMP/jeter.build
then
	wget -O - "https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz" |\\
		gunzip -c |\\
		jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP convertliftoverchain -R2 ${dict} > TMP/jeter.chain

	bigBedToBed -udcDir=TMP/CACHE TMP/jeter.bb stdout |\\
		jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP bedliftover --chain TMP/jeter.chain -R ${dict} |\\
        	LC_ALL=C sort -t '\t' -k1,1 -k2,2n -T TMP |\\
	        bgzip > TMP/${TAG}.bed.gz

else
	bigBedToBed -udcDir=TMP/CACHE TMP/jeter.bb stdout |\\
	        jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP bedrenamechr -R ${fasta} --column 1 --convert SKIP  |\\
        	LC_ALL=C sort -t '\t' -k1,1 -k2,2n -T TMP |\\
	        bgzip > TMP/${TAG}.bed.gz
fi


tabix -p bed -f TMP/${TAG}.bed.gz

mv TMP/${TAG}.bed.gz ./
mv TMP/${TAG}.bed.gz.tbi ./

echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="${whatis}  jvarkit.info.bean=(name|score)">' > ${TAG}.header

cat << EOF > ${TAG}.md
Heart enhancer predictions (Dickel et al.) ${url}. These predicted, putative heart enhancers are inferred by combining the ChIP-seq enrichment for either H3K27ac or p300 across 1 human fetal sample [1], 
8 human postnatal [1,2], 10 mouse fetal [3-7] and 9 mouse postnatal [3,5-8]). 
EOF
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
    def TAG = "ENHANCER_DICKEL"
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
"""
}
