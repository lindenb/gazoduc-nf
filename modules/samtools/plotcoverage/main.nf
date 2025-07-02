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

include {moduleLoad;getVersionCmd} from '../../modules/utils/functions.nf'

process SAMTOOLS_DEPTH_PLOT_COVERAGE {
tag "${row.bam}"
afterScript "rm -rf TMP"
cpus 1
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta4),path(gtf),path(gtfidx)
	tuple val(meta ),path("BAMS/*")
output:
	tuple val(meta),val("*.pdf"),emit:pdf
	path("versions.yml"),emit:versions
script:
	def prefix = row.prefix?:""
	def contig = row.contig
	def start = row.start
	def end = row.end
	def mapq = row.mapq?:1

"""
hostname 1>&2
set -x
mkdir -p TMP
find BAMS/ -name "*.bam" -o -name "*.cram" | sort > TMP/bams.list


echo "SN\tDEPTH\tmaxDP" > TMP/depths.tsv

# call samtools depth
cat "TMP/bams.list" | samtools samples  | cut -f1,2 | while read SN BAM
do
	
	samtools depth  -a -r "${contig}:${start}-${end}" "\${BAM}" | cut -f 2,3 > TMP/\${SN}.depth.txt
	echo -n "\${SN}\tTMP/\${SN}.depth.txt\t" >> TMP/depths.tsv
	# get max cov
	cut -f2 TMP/\${SN}.depth.txt | sort -T . -n | tail -1 >> TMP/depths.tsv
done

# exons
echo "contig\tstart\tend" > TMP/exons.bed
if ${meta.containsKey("gtf")} ; then
	tabix "${meta.gtf}"  "${contig}:${start}-${end}" |\
	awk -F '\t' '(\$3=="exon")' |\
	cut -f1,4,5 | sort -T TMP | uniq >> TMP/exons.bed
fi


cat "${moduleDir}/plot.R" |\\
	m4 -D__CONTIG__=${contig} \\
		-D__START__=${start} \\
		-D__END__=${end}  > TMP/jeter.R

R --vanilla < TMP/jeter.R


"""
}
