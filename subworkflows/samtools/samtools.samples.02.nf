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



include {SQRT_FILE} from '../../modules/utils/sqrt.nf'
include {assertFileExists;isBlank;moduleLoad;parseBoolean} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {ALL_REFERENCES} from '../../modules/reference/all.references.01.nf'

workflow SAMTOOLS_SAMPLES02 {
	take:
		meta
		genomeId /* or empty string */
		bams
	main:

		for(String s : ["allow_multiple_references","with_header","allow_duplicate_samples"]) {
			if(!meta.containsKey(s)) throw new IllegalArgumentException("meta."+s+" undefined");
			}

		version_ch = Channel.empty()
		
		sqrt_ch = SQRT_FILE([min_file_split:100, suffix:".list"],bams)
		version_ch = version_ch.mix(sqrt_ch.version)

		all_refs_ch = ALL_REFERENCES()
		version_ch = version_ch.mix(all_refs_ch.version)

		chunk_ch = sqrt_ch.output.splitText().map{it.trim()}

		sn_ch = ST_SAMPLE(all_refs_ch.output,chunk_ch)
		version_ch = version_ch.mix(sn_ch.version)
	
		digest_ch = DIGEST(meta, genomeId, sn_ch.output.collect())
		version_ch = version_ch.mix(digest_ch.version)

		version_ch = MERGE_VERSION("samtools samples", version_ch.collect())		
	emit:
		references = all_refs_ch.output
		version = version_ch
		output = digest_ch.output
	}


process ST_SAMPLE {
tag "${bams.name}"
afterScript "rm -rf TMP"
input:
	path(references)
	path(bams)
output:
	path("sample2bam.tsv"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 2>&1
${moduleLoad("samtools")}
set -o pipefail


samtools samples  -F "${references}" < "${bams}" > sample2bam.tsv

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">extract sample names from BAM metadata</entry>
	<entry key="input">${bams}</entry>
        <entry key="samtools.version">\$(samtools  --version | head -n 1| cut -d ' ' -f2)</entry>
</properties>
EOF
"""
}

process DIGEST {
executor "local"
tag "N=${L.size()}"
input:
	val(meta)
	val(gemomeId)
	val(L)
output:
	path("sample2bam.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def fasta = isBlank(genomeId) || genomeId.equals(".") ? "" : params.genomes[genomeId].fasta
	def allow_multiple_references = (meta.allow_multiple_references as boolean)
"""
hostname 1>&2
set -o pipefail

cat ${L.join( " ") } |\
	${allow_multiple_references || isBlank(fasta)?"":" awk -F '\t' '(\$3 == "${fasta}")' |"} \
	sort -T . -t '\t' -k1,1  |\
	uniq |\
	awk -F '\t' 'function uniq(S) {P=S;i=1; while(P in U) {i++;P=sprintf("%s.%d",S,i);} U[P]++; return P;} {printf("%s\t%s\t%s\t%s\\n",\$1,uniq(\$1),\$2,\$3);}' > jeter.tsv

# no empty samples
awk -F '\t' '(\$1==".")' jeter.tsv > jeter.txt 
test ! -s jeter.txt

# no empty ref
awk -F '\t' '(\$4==".")' jeter.tsv > jeter.txt 
test ! -s jeter.txt

# no dup samples
if  ${!parseBoolean(meta.allow_duplicate_samples)}  ; then
	cut -f 1 jeter.tsv | sort -T . | uniq -d > jeter.txt
	test ! -s jeter.txt
fi

# no mulitple refs
if  ${!allow_multiple_references}  ; then
	test `cut -f 4 jeter.tsv | sort -T . | uniq  | wc -l ` == 1
fi


# no dup bam
cut -f 3 jeter.tsv | sort -T . | uniq -d > jeter.txt
test ! -s jeter.txt

awk -F '\t'  '(\$1!=\$2)' jeter.tsv > duplicate.names.txt


if  ${parseBoolean(meta.with_header)} ; then
	echo "sample\tnew_sample\tbam\treference" > sample2bam.tsv
fi

cat jeter.tsv >>  sample2bam.tsv


##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">extract sample names from BAM metadata</entry>
	<entry key="input.count">${L.size()}</entry>
        <entry key="n-samples">\$(wc -l < sample2bam.tsv )</entry>
        <entry key="samples">\$(cut -f 1 sample2bam.tsv |paste -s -d ' ')</entry>
        <entry key="duplicate names"><pre>\$(cut -f1,3 duplicate.names.txt)</pre></entry>
        <entry key="distinct references"><pre>\$(cut -f4 sample2bam.tsv | sort | uniq)</pre></entry>
</properties>
EOF
"""
}
