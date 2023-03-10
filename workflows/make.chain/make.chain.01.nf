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
nextflow.enable.dsl=2

params.source_fasta  = "NO_FILE"
params.dest_fasta  = "NO_FILE"

include {REPEAT_MASKER_FASTA as MASK1; REPEAT_MASKER_FASTA as MASK2} from '../../subworkflows/repeat.masker/repeat.mask.fasta.nf'

//include {SAMTOOLS_SAMPLES_01} from '../../subworkflows/samtools/samtools.samples.01.nf'
//include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {moduleLoad;runOnComplete} from '../../modules/utils/functions.nf'
//include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
//include {SIMPLE_ZIP_01} from '../../modules/utils/zip.simple.01.nf'



workflow {
	ch1 = MAKE_CHAIN(params, params.source_fasta, params.dest_fasta )
	}

runOnComplete(workflow)

workflow MAKE_CHAIN {
	take:
		meta
		source_fasta
		dest_fasta
	main:
		version_ch  = Channel.empty()
		src_ch = MASK1([extraRepeatMasker:""],source_fasta)
		dest_ch = MASK2([extraRepeatMasker:""],dest_fasta)

		fa2bits_ch = FA2BITS([:],src_ch.fasta.map{T->["src",T]}.mix(dest_ch.fasta.map{T->["dest",T]}) )
		
		src2bit_ch = fa2bits_ch.output.filter{T->T[0].equals("src")}.map{T->T[1]}.first()
		src2tab_ch = fa2bits_ch.output.filter{T->T[0].equals("src")}.map{T->T[2]}.first()
		dest2bit_ch = fa2bits_ch.output.filter{T->T[0].equals("dest")}.map{T->T[1]}.first()
		dest2tab_ch = fa2bits_ch.output.filter{T->T[0].equals("dest")}.map{T->T[2]}.first()

		splitnew_ch = SPLIT_NEW([:], dest_ch.fasta)
		
		align_ch = ALIGN([:],  src2bit_ch, dest2bit_ch , splitnew_ch.lft, splitnew_ch.output.splitText().map{it.trim()} )
		
		join_ch = CHAIN_MERGE_SORT([:], align_ch.output.collect())
		
	}


process FA2BITS {
tag "${fasta.name} ${type}"
input:
	val(meta)
	tuple val(type),path(fasta)
output:
	tuple val(type),path("${fasta.simpleName}.2bit"),path("${fasta.simpleName}.tab"),emit:output
script:
"""
hostname 1>&2
${moduleLoad("ucsc")}

faToTwoBit "${fasta}" "${fasta.simpleName}.2bit"
twoBitInfo "${fasta.simpleName}.2bit" "${fasta.simpleName}.tab"
"""
}


process SPLIT_NEW {
tag "${fasta.name}"
input:
	val(meta)
	path(fasta)
output:
	path("subseq.list"),emit:output
	path("output.lft"),emit:lft
script:
	def splitSize = 5000
"""
hostname 1>&2
${moduleLoad("ucsc")}

mkdir -p OUT

faSplit -prefixLength=5 -verbose=3 -outDirDepth=3 -lift=output.lft  size "${fasta}" ${splitSize}  OUT/
find \${PWD}/OUT -type f -name "*.fa" > subseq.list
test -s subseq.list
"""
}

process ALIGN {
tag "${src2bit.name} -> ${dest2bit.name} (${file(newfa).name})"
afterScript "rm -rf TMP"
input:
	val(meta)
	path(src2bit)
	path(dest2bit)
	path(lft)
	val(newfa)
output:
	path("output.chain"),emit:output
script:
	//def extraBlat = "-tileSize=12 -minScore=100 -minIdentity=98 -fastMap"
	def extraBlat = ""
"""
hostname 1>&2
${moduleLoad("ucsc")}

mkdir -p TMP

lastz "${src2bit}" "${newfa}" > TMP/output.lav

lavToPsl TMP/output.lav TMP/lastz.psl

liftUp -pslQ TMP/output.psl "${lft}" warn TMP/lastz.psl

axtChain -linearGap=medium -psl TMP/output.psl "${src2bit}" "${dest2bit}" output.chain

"""
}

process CHAIN_MERGE_SORT {
tag "N=${L.size()}"
input:
	val(meta)
	val(L)
output:
	path("chain"),emit:dir
script:
"""
hostname 1>&2
${moduleLoad("ucsc")}

chainMergeSort ${L.join(" ")} |\
	chainSplit chain stdin
"""
}

