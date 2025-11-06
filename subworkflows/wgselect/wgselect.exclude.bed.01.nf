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

include { SCATTER_TO_BED } from '../../subworkflows/picard/picard.scatter2bed.nf'
include {k1_signature} from '../../modules/utils/k1.nf'

def k1 = k1_signature()

workflow WGSELECT_EXCLUDE_BED_01 {
	take:
		genome
	main:
		to_merge_ch = Channel.empty()

		gaps_ch = SCATTER_TO_BED(genome)
		to_merge_ch = to_merge_ch.mix(gaps_ch.output)

		if(params.wgselect.with_rmsk.toBoolean()) {
			rmsk_ch = RMSK(genome)
			to_merge_ch = to_merge_ch.mix(rmsk_ch.bed)
			}

		if(params.wgselect.with_encode_exclude.toBoolean()) {
			x2_ch = EXCLUDE_ENCODE(genome)
			to_merge_ch = to_merge_ch.mix(x2_ch.bed)
			}


		if(params.wgselect.with_lcr.toBoolean()) {
			x3_ch = LOW_COMPLEXITY_REGIONS(genome)
			to_merge_ch = to_merge_ch.mix(x3_ch.bed)
			}

		if(params.wgselect.with_simple_repeats .toBoolean()) {
			x4_ch = SIMPLE_REPEATS(genome)
			to_merge_ch = to_merge_ch.mix(x4_ch.bed)
			}
		all_x_ch = MERGE_REGIONS(to_merge_ch.collect())
	emit:
		bed = all_x_ch.bed
	}




process RMSK {
label "process_short"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	path(genome)
output:
	path("rmsk.bed"), emit:bed
script:
	def fai = genome.find{it.name.endsWith(".fai")}
	def fasta = genome.find{it.name.endsWith("a")}
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP

cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter.a
1:${k1.hg19}\thttps://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz
1:${k1.hg38}\thttps://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
EOF

awk '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sed 's/^chr//' |  sort -T TMP -t '\t' -k1,1 > TMP/jeter.b

URL=`join -t '\t' -1 1 -2 1 -o 1.2 TMP/jeter.a TMP/jeter.b`


wget -O - "\${URL}" |\\
	gunzip -c |\\
	cut -f6-8 |\\
	jvarkit bedrenamechr -f "${fasta}" --column 1 --convert SKIP |\
		LC_ALL=C sort -t '\t' -T . -k1,1 -k2,2n > rmsk.bed

test -s rmsk.bed
"""
}

process EXCLUDE_ENCODE {
label "process_short"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	path(genome)
output:
	path("excude.encode.bed"),emit:bed
script:
        def fai = genome.find{it.name.endsWith(".fai")}
        def fasta = genome.find{it.name.endsWith("a")}

"""
hostname 1>&2
set -o pipefail
mkdir -p TMP

cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter.a
1:${k1.hg19}\thttps://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz
1:${k1.hg38}\thttps://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz?raw=true
EOF

awk '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sed 's/^chr//' |  sort -T TMP -t '\t' -k1,1 > TMP/jeter.b

URL=`join -t '\t' -1 1 -2 1 -o 1.2 TMP/jeter.a TMP/jeter.b`


wget -O - "\${URL}" |\
	gunzip -c |\
	cut -f1-3 |\
	jvarkit bedrenamechr -f "${fasta}" --column 1 --convert SKIP |\
	LC_ALL=C sort -t '\t' -T . -k1,1 -k2,2n > excude.encode.bed

test -s excude.encode.bed
"""
}


process LOW_COMPLEXITY_REGIONS {
label "process_short"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	path(genome)
output:
	path("lcr.bed"), emit:bed
script:
	def fai = genome.find{it.name.endsWith(".fai")}
	def fasta = genome.find{it.name.endsWith("a")}
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP

cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter.a
1:${k1.hg19}\thttps://github.com/lh3/varcmp/blob/master/scripts/LCR-hs37d5.bed.gz?raw=true
1:${k1.hg38}\thttps://github.com/lh3/varcmp/blob/master/scripts/LCR-hs38.bed.gz?raw=true
EOF

awk '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sed 's/^chr//' |  sort -T TMP -t '\t' -k1,1 > TMP/jeter.b

URL=`join -t '\t' -1 1 -2 1 -o 1.2 TMP/jeter.a TMP/jeter.b`
		
wget -O - "\${URL}" |\
	gunzip -c |\
	cut -f1-3 |\
	jvarkit bedrenamechr -f "${fasta}" --column 1 --convert SKIP |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n > lcr.bed

test -s lcr.bed
"""
}

process SIMPLE_REPEATS {
label "process_short"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	path(genome)
output:
	path("simple_repeats.bed"), emit:bed
script:
        def fai = genome.find{it.name.endsWith(".fai")}
        def fasta = genome.find{it.name.endsWith("a")}

"""
hostname 1>&2
set -o pipefail
mkdir -p TMP

cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter.a
1:${k1.hg19}\thttps://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz
1:${k1.hg38}\thttps://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz
EOF

awk '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sed 's/^chr//' |  sort -T TMP -t '\t' -k1,1 > TMP/jeter.b

URL=`join -t '\t' -1 1 -2 1 -o 1.2 TMP/jeter.a TMP/jeter.b`


wget -O - "\${URL}" |\
	gunzip -c |\
	cut -f2-4 |\
	jvarkit bedrenamechr -f "${fasta}" --column 1 --convert SKIP |\
		LC_ALL=C sort -t '\t' -T . -k1,1 -k2,2n > simple_repeats.bed


test -s simple_repeats.bed
"""
}


process MERGE_REGIONS {
label "process_short"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	path("BED/*")
output:
	path("bad_regions.bed"),emit:bed
script:
"""
hostname 1>&2
set -o pipefail
cut -f1,2,3  BED/*.bed |\\
LC_ALL=C sort  -T . -k1,1 -k2,2n  |\\
	bedtools merge > bad_regions.bed
"""
}


