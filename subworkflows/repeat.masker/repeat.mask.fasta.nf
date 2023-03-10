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
include {moduleLoad;getVersionCmd} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'


/**
 * APPLY RepeatMasker to a reference genome
 *
 */
workflow REPEAT_MASKER_FASTA {
	take:
		meta
		fasta /* indexed with samtools */
	main:
		version_ch = Channel.empty()
		chroms_ch = INIT_AND_LIST_CONTIGS(meta, fasta)
		version_ch = version_ch.mix(chroms_ch.version)

		each_contig = chroms_ch.output.splitText().map{it.trim()}

		mask_ch  = MASK_FASTA(meta, fasta, each_contig)
		version_ch = version_ch.mix(mask_ch.version)

		join_ch = JOIN(meta, fasta, chroms_ch.output, mask_ch.output.map{T->T.join("\t")}.collect())
		version_ch = version_ch.mix(join_ch.version)

		version_ch = MERGE_VERSION(meta, "repeatmasker", "repeat masker fasta", version_ch.collect())
	emit:
		fasta = join_ch.output
		version = version_ch
	}

process INIT_AND_LIST_CONTIGS {
tag "${file(reference).name}"
executor "local"
maxForks 1 /* we'll fill RepatMasker Cache for the first time here */
input:
	val(meta)
	val(reference)
output:
	path("contigs.txt"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("repeatmasker")}
set -x
mkdir -p TMP
cut -f 1 "${reference}.fai" > contigs.txt
test -s contigs.txt


#
# on our cluster test this directory is a symbolic line that goes to elswhere because I've got a small quotas on the login node 
#
test -L "\${HOME}/.RepeatMaskerCache"

# better than nothing... try to create a LOCK file
test ! -f "\${HOME}/.RepeatMaskerCache/LOCK"
touch "\${HOME}/.RepeatMaskerCache/LOCK"

#
# run repeat masker once to fill the cache and avoid collisition later when running each contig in parallel
#

echo ">TMP" > TMP/jeter.fa
tr "\\0" "A" < /dev/zero | fold -w 60 | head -n 100 >> TMP/jeter.fa


RepeatMasker -dir TMP  ${meta.extraRepeatMasker?:""} TMP/jeter.fa

rm "\${HOME}/.RepeatMaskerCache/LOCK"

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">extract contigs from reference</entry>
	<entry key="reference">${reference}</entry>
	<entry key="repeat.masker.version">\$(RepeatMasker -v)</entry>
</properties>
EOF
"""
}

process MASK_FASTA {
tag "${file(reference).name} ${contig}"
cpus 4
afterScript "rm -rf TMP"

input:
	val(meta)
	val(reference)
	val(contig)
output:
	tuple val(contig),path("${contig}.fa.masked.gz"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("samtools repeatmasker")}
set -o pipefail


mkdir -p TMP

samtools faidx "${reference}" "${contig}" > TMP/${contig}.fa
RepeatMasker  -parallel ${task.cpus} -dir TMP ${meta.extraRepeatMasker?:""} TMP/${contig}.fa 1>&2
gzip --best  TMP/${contig}.fa.masked
mv TMP/${contig}.fa.masked.gz ./

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">apply repeat masker to contig</entry>
	<entry key="repeat.masker.version">\$(RepeatMasker -v)</entry>
	<entry key="reference">${reference}</entry>
	<entry key="contig">${contig}</entry>
	<entry key="version">${getVersionCmd("samtools")}</entry>
</properties>
EOF
"""
}


process JOIN {
tag "N=${L.size()}"
input:
	val(meta)
	val(reference)
	val(contigs)
	val(L)

output:
	path("${file(reference).getBaseName()}.masked.fa"),emit:output
	path("${file(reference).getBaseName()}.masked.fa.fai")
	path("${file(reference).getBaseName()}.masked.dict")
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("samtools")}
mkdir -p TMP

awk '{printf("%s\t%d\\n",\$1,NR);}' '${contigs}' | sort -t '\t' -T TMP -k1,1 > TMP/jeter.a

cat << EOF | sort -t '\t' -T TMP -k1,1 > TMP/jeter.b
${L.join("\n")}
EOF

join -t '\t' -1 1 -2 1 TMP/jeter.a TMP/jeter.b | sort -t '\t' -T TMP -k2,2n | cut -f 3 > TMP/jeter.list

xargs -a TMP/jeter.list -L 10 cat | gunzip -c > TMP/jeter.fa

mv TMP/jeter.fa "${file(reference).getBaseName()}.masked.fa"

samtools faidx "${file(reference).getBaseName()}.masked.fa"

samtools dict --alias -o "${file(reference).getBaseName()}.masked.dict"  "${file(reference).getBaseName()}.masked.fa" 


#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">merge masked contigs</entry>
	<entry key="version">${getVersionCmd("samtools awk")}</entry>
</properties>
EOF
"""
}
