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

process SAMTOOLS_DEPTH_PLOT_COVERAGE {
tag "${meta.id?:""}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(gtf),path(gtfidx)
	tuple val(meta),path(fasta),path(fai) ,path(dict),path("BAMS/*")
output:
	tuple val(meta),path("*.pdf"),emit:pdf
	tuple val(meta),path("*.html"),emit:html
	path("versions.yml"),emit:versions
script:
	def has_gtf =gtf?true:false
    def contig = meta.contig
    def start = (meta.start as int)
    def end = (meta.end as int)
    def factor = ((task.ext.factor?:5) as double)
    def len = (end-start+1)
    def mid = start + len/2;
    def xlen = ((len*factor) as int)
    def xstart = java.lang.Math.max(1,((mid-xlen) as int))
    def xend = mid+xlen
    def title = task.ext.title?:meta.title
	def mapq = task.ext.mapq?:20
	def prefix = task.ext.prefix?:meta.id
"""
hostname 1>&2
set -x
mkdir -p TMP
find BAMS/ -name "*.bam" -o -name "*.cram" | sort > TMP/bams.list


echo "SN\tDEPTH\tmaxDP" > TMP/depths.tsv

# call samtools depth
cat "TMP/bams.list" | samtools samples  | cut -f1,2 | sort -t '\t' -k1,1V | while read SN BAM
do
	samtools depth --threads ${task.cpus} --min-MQ ${mapq} -a -r "${contig}:${xstart}-${xend}" "\${BAM}" |\\
		cut -f 2,3 > TMP/\${SN}.depth.txt
	echo -n "\${SN}\tTMP/\${SN}.depth.txt\t" >> TMP/depths.tsv
	# get max cov
	cut -f2 TMP/\${SN}.depth.txt | LC_ALL=C sort -T . -n | tail -1 >> TMP/depths.tsv
done

# exons
echo "contig\tstart\tend" > TMP/exons.bed
if ${has_gtf}
then
	tabix "${gtf}"  "${contig}:${xstart}-${xend}" |\\
		awk -F '\t' '(\$3=="exon")' |\\
		cut -f1,4,5 |\\
		LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\\
		bedtools merge  >> TMP/exons.bed
fi


cat "${moduleDir}/plot.R" |\\
	m4 -D__CONTIG__=${contig} \\
		-D__START__=${start} \\
		-D__END__=${end} \\
		-D__LEN__=${len} \\
		-D__TITLE__="${title}" \\
		-D__PREFIX__="${prefix}" \\
		-D__XSTART__=${xstart} \\
		-D__XEND__=${xend} \\
		-D__XLEN__=${xlen} > TMP/jeter.R

R --vanilla < TMP/jeter.R

cat << EOF > TMP/jeter.html
<tr>
	<td><code>${contig}:${start}-${end}</td>
	<td>${len}</td>
	<td>${title}</td>
	<td><a href="${prefix}.pdf">${prefix}.pdf</a></td>
</tr>
EOF

mv TMP/jeter.html ${prefix}.html

cat << END_VERSIONS > versions.yml
"${task.process}":
	samtools: "\$(samtools version | awk '(NR==1) {print \$NF;}')"
	R: todo
END_VERSIONS
"""
}
