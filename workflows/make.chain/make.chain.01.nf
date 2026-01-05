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
nextflow.enable.dsl=2


def gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults()


gazoduc.build("source_fasta","NO_FILE").
	desc("OLD fasta genome. Must be indexed wit faidx").
	required().
	existingFile().
	put()

gazoduc.build("source_is_repeat_masked",false).
	desc("true if source is hard repeat masked").
	setBoolean().
	put()


gazoduc.build("dest_fasta","NO_FILE").
	desc("NEW fasta genome. Must be indexed wit faidx").
	required().
	existingFile().
	put()

gazoduc.build("dest_is_repeat_masked",false).
	desc("true if destination is hard repeat masked").
	setBoolean().
	put()


gazoduc.build("split_size",10000).
        desc("Split the genome into fragments of 'x' bp.").
	setInteger().
        put()

gazoduc.build("n_submit",10000).
        desc("Submit a batch 'x' sequences when running blat/lastz.").
	setInteger().
        put()


include {REPEAT_MASKER_FASTA as MASK1; REPEAT_MASKER_FASTA as MASK2} from '../../subworkflows/repeat.masker/repeat.mask.fasta.nf'
//include {SAMTOOLS_SAMPLES_01} from '../../subworkflows/samtools/samtools.samples.01.nf'
//include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {moduleLoad;runOnComplete;parseBoolean} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
//include {SIMPLE_ZIP_01} from '../../modules/utils/zip.simple.01.nf'


if( params.help ) {
    gazoduc.usage().
	name("make liftover chain").
	desc("make liftover chains from two fasta sequences").
	print();
    exit 0
} else {
   gazoduc.validate();
}



workflow {
	ch1 = MAKE_CHAIN(params, file(params.source_fasta), file(params.dest_fasta) )
	}

runOnComplete(workflow)

workflow MAKE_CHAIN {
	take:
		meta
		source_fasta
		dest_fasta
	main:
		version_ch  = Channel.empty()

		if(parseBoolean(meta.source_is_repeat_masked)) {
			masked_source = Channel.fromPath(source_fasta)
			}
		else
			{
			src_ch = MASK1([extraRepeatMasker:""],source_fasta)
			version_ch = version_ch.mix(src_ch.version)
			masked_source = dest_ch.fasta
			}

		if(parseBoolean(meta.dest_is_repeat_masked)) {
			masked_dest = Channel.fromPath(dest_fasta)
			}
		else
			{
			dest_ch = MASK2([extraRepeatMasker:""],dest_fasta)
			version_ch = version_ch.mix(dest_ch.version)
			masked_dest = dest_ch.fasta
			}
	

		fa2bits_ch = FA2BITS([:],masked_source.map{T->["src",T]}.mix(masked_dest.map{T->["dest",T]}) )
		version_ch = version_ch.mix(fa2bits_ch.version)	

		src2bit_ch = fa2bits_ch.output.filter{T->T[0].equals("src")}.map{T->T[1]}.first()
		src2tab_ch = fa2bits_ch.output.filter{T->T[0].equals("src")}.map{T->T[2]}.first()
		dest2bit_ch = fa2bits_ch.output.filter{T->T[0].equals("dest")}.map{T->T[1]}.first()
		dest2tab_ch = fa2bits_ch.output.filter{T->T[0].equals("dest")}.map{T->T[2]}.first()

		splitnew_ch = SPLIT_NEW([:], dest_ch.fasta)
		version_ch = version_ch.mix(splitnew_ch.version)		

		align_ch = ALIGN([:],  src2bit_ch, dest2bit_ch , splitnew_ch.lft, splitnew_ch.output.splitText().map{it.trim()} )
		version_ch = version_ch.mix(align_ch.version)

		join_ch = CHAIN_MERGE_SORT(
				meta,
				source_fasta, dest_fasta,
				src2tab_ch, dest2tab_ch,
				align_ch.output.collect()
				)
		version_ch = version_ch.mix(join_ch.version)

		version_ch = MERGE_VERSION(meta, "makechain", "make chain",version_ch.collect())

	emit:
			output = join_ch.output
			version = version_ch		
	}


process FA2BITS {
tag "${fasta.name} ${type}"
input:
	val(meta)
	tuple val(type),path(fasta)
output:
	tuple val(type),path("${fasta.simpleName}.2bit"),path("${fasta.simpleName}.tab"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("ucsc")}

faToTwoBit "${fasta}" "${fasta.simpleName}.2bit"
twoBitInfo "${fasta.simpleName}.2bit" "${fasta.simpleName}.tab"

cat <<- EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">convert fasta to two bits format</entry>
</properties>
EOF
"""
}



process SPLIT_NEW {
tag "${fasta.name}"
afterScript "rm -rf TMP"
input:
	val(meta)
	path(fasta)
output:
	path("subseq.list"),emit:output
	path("output.lft"),emit:lft
	path("version.xml"),emit:version
script:
	def splitSize = params.split_size?:1000000
	def nLines = params.n_submit?:1000
	def maxN = ( (splitSize as int)*0.8) as int)
"""
hostname 1>&2
${moduleLoad("ucsc")}

mkdir -p TMP
mkdir -p OUT

faSplit -prefixLength=5  -maxN=${maxN} -verbose=3 -outDirDepth=3 -lift=output.lft  size "${fasta}" ${splitSize}  OUT/
find \${PWD}/OUT -type f -name "*.fa" | split -a 8 --lines=${nLines} --additional-suffix=.list - OUT/cluster > subseq.list

find \${PWD}/OUT -type f -name "cluster*.list" > subseq.list
test -s subseq.list

cat <<- EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">split reference into smaller parts</entry>
</properties>
EOF
"""
}

process ALIGN {
tag "${src2bit.name} -> ${dest2bit.name} (${file(fastalist).name})"
afterScript "rm -rf TMP"
memory "10g"
input:
	val(meta)
	path(src2bit)
	path(dest2bit)
	path(lft)
	val(fastalist)
output:
	path("output.chain"),emit:output
	path("version.xml"),emit:version
script:
	//def extraBlat = "-tileSize=12 -minScore=100 -minIdentity=98 -fastMap"
	def extraBlat = ""
"""
hostname 1>&2
${moduleLoad("ucsc")}

mkdir -p TMP

xargs -a "${fastalist}" -L 10 cat > TMP/jeter.fa

lastz "${src2bit}" "TMP/jeter.fa" > TMP/output.lav

lavToPsl TMP/output.lav TMP/lastz.psl

liftUp -pslQ TMP/output.psl "${lft}" warn TMP/lastz.psl

axtChain -linearGap=medium -psl TMP/output.psl "${src2bit}" "${dest2bit}" output.chain


cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">align on genome</entry>
</properties>
EOF
"""
}

process CHAIN_MERGE_SORT {
//tag "N=${L.size()}"
afterScript "rm -rf TMP"
input:
        val(meta)
        path(src_fa)
        path(dest_fa)
        path(src_sizes)
        path(dest_sizes)
        val(L)
output:
	path("${meta.prefix?:""}${src_fa.simpleName}_To_${dest_fa.simpleName}.chain"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2

${moduleLoad("ucsc")}


chainMergeSort ${L.join(" ")} | chainSplit TMP stdin

cat TMP/*.chain > TMP/jeter.txt

chainSort TMP/jeter.txt TMP/all.sorted.chain

chainNet TMP/all.sorted.chain "${src_sizes}" "${dest_sizes}" TMP/output.net /dev/null

netChainSubset TMP/output.net TMP/all.sorted.chain "${meta.prefix?:""}${src_fa.simpleName}_To_${dest_fa.simpleName}.chain"

cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">merge final chain</entry>
</properties>
EOF
"""
}

