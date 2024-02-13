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



include {SQRT_FILE} from '../../modules/utils/sqrt.nf'
include {assertFileExists;isBlank;moduleLoad;parseBoolean} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {ALL_GENOMES} from '../../modules/reference/all.genomes.01.nf'

workflow SAMTOOLS_SAMPLES {
	take:
		bams
	main:

		version_ch = Channel.empty()
		
		sqrt_ch = SQRT_FILE([min_file_split:100, suffix:".list"],bams)
		version_ch = version_ch.mix(sqrt_ch.version)

		all_genomes_ch = ALL_GENOMES()
		version_ch = version_ch.mix(all_genomes_ch.version)

		all_genomes_rows_ch = all_genomes_ch.output.splitCsv(sep:'\t',header:true)

		all_refs_file = all_genomes_rows_ch.
			map{T->T.fasta}.
			collectFile(name: 'references.txt', newLine: true)

		chunk_ch = sqrt_ch.output.splitText().map{it.trim()}.combine(all_refs_file)

		sn_ch = ST_SAMPLE(chunk_ch)
		version_ch = version_ch.mix(sn_ch.version)
	
		digest_ch = DIGEST(sn_ch.output.collect())
		version_ch = version_ch.mix(digest_ch.version)

		rows_ch = digest_ch.output.splitCsv(sep:'\t',header:true).
				combine(all_genomes_rows_ch).
				filter{T->T[0].fasta.equals(T[1].fasta)}.
				map{T->T[0].plus(T[1])}.
				map{T->{
					if(T.indexed.equals("T")) T.put("indexed":true);
					else if(T.indexed.equals("F")) T.put("indexed":false);
					return T;
					}}.
				map{T->T.plus("reference":T.fasta)}

		version_ch = MERGE_VERSION("samtools samples", version_ch.collect())		
	emit:
		genomes = all_genomes_ch.output
		rows = rows_ch
		version = version_ch
		output = digest_ch.output
	}


process ST_SAMPLE {
tag "${file(bams).name}"
afterScript "rm -rf TMP"
input:
	tuple val(bams),val(references)
output:
	path("sample2bam.tsv"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 2>&1
${moduleLoad("samtools")}
set -o pipefail

samtools samples  -i -F "${references}" < "${bams}" > sample2bam.tsv

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
afterScript "rm -rf TMP"
input:
	val(L)
output:
	path("sample2bam.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def allow_duplicate_samples = task.ext.allow_duplicate_samples
"""
hostname 1>&2
set -o pipefail
set -x

mkdir -p TMP

cat ${L.join( " ") } | sort -T TMP | uniq > TMP/jeter2.tsv

sort -T TMP -t '\t' -k1,1  TMP/jeter2.tsv |\
	uniq |\
	awk -F '\t' 'function uniq(S) {P=S;i=1; while(P in U) {i++;P=sprintf("%s.%d",S,i);} U[P]++; return P;} {printf("%s\t%s\t%s\t%s\t%s\\n",\$1,uniq(\$1),\$2,\$3,\$4);}' > TMP/jeter.tsv
rm TMP/jeter2.tsv


# no dup samples
if  ${!allow_duplicate_samples}  ; then
	cut -f 1 TMP/jeter.tsv | sort -T TMP | uniq -d > TMP/jeter.txt
	cat TMP/jeter.txt 1>&2
	test ! -s TMP/jeter.txt
fi


# no dup bam
cut -f 3 TMP/jeter.tsv | sort -T TMP | uniq -d > TMP/jeter.txt
cat TMP/jeter.txt 1>&2
test ! -s TMP/jeter.txt

awk -F '\t'  '(\$1!=\$2)' TMP/jeter.tsv > duplicate.names.txt


echo "sample\tnew_sample\tbam\tindexed\tfasta" > sample2bam.tsv

cat TMP/jeter.tsv >>  sample2bam.tsv


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
