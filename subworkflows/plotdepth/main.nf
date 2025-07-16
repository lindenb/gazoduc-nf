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
include {SAMTOOLS_SAMPLES02} from '../../subworkflows/samtools/samtools.samples.02.nf'


workflow PLOT_COVERAGE_01 {
	take:
		meta
		fasta
		fai
		dict
		gtf
		bed
		bams //[bam,bai]
	main:
		versions = Channel.empty()
		EXTEND_BED(fasta,fai,dict, bed)
		versions = versions.mix(EXTEND_BED.out.versions)

		c1_ch = EXTEND_BED.out.bed.
			map{it[1]}.
			splitCsv(header: true,sep:'\t',strip:true)



		DRAW_COVERAGE(fasta,fai,dict,gtf,c1_ch.combine(bams))
		versions = versions.mix(DRAW_COVERAGE.out.versions)

		MERGE_PDFS(
			DRAW_COVERAGE.out.pdf
				.map{T->[[id:T[0].chrom+"_"+T[0].delstart+"_"+T[0].delend],T[1]]}
				.groupTuple()
			)
		versions = versions.mix(MERGE_PDFS.out.versions)


		ZIP_ALL(MERGE_PDFS.out.pdf.map{it[1]}.collect().map{[meta,it]})
		versions = versions.mix(ZIP_ALL.out.versions)
	emit:
		zip =ZIP_ALL.out.zip
		versions
	}


process EXTEND_BED {
	label "process_single"
	tag "${meta.id} ${bed.name}"
	conda "${moduleDir}/../../conda/bioinfo.01.yml"
	input:
		tuple val(meta1),path(fasta)
		tuple val(meta2),path(fai)
		tuple val(meta3),path(dict)
		tuple val(meta),path(bed)
	output:
		tuple val(meta),path("*.bed"),emit:bed
		path("versions.yml"),emit:versions
	script:
		def extend = task.ext.extend_bed?:"3.0"
		def max_sv_len = ((task.ext.max_sv_length?:-1) as int)
		def min_sv_len = ((task.ext.min_sv_length?:-1) as int)
	"""
	hostname 1>&2
	mkdir TMP
	cut -f1,2 "${fai}" > TMP/jeter.genome

cat << 'EOF' > jeter.awk
BEGIN {
	M=${max_sv_len};
	N=${min_sv_len};
	}

/^(browser|track|#)/ {
	next;
	}

	{
	L=int(\$3)-int(\$2);
	if(M>0 && L>M) next;
	if(N>0 && L<N) next;
	printf("%s\t%s\t%s\t%s\\n",\$1,\$2,\$3, (NF==3 || \$4==""?".":\$4) );
	}
EOF

	jvarkit  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP bedrenamechr -R "${fasta}" -c 1 "${bed}" | \\
		awk -F '\t' -f jeter.awk > TMP/jeter1.bed

	cut -f 1,2,3 TMP/jeter1.bed |\\
		bedtools slop -i - -g TMP/jeter.genome ${extend.toString().contains(".")?"-pct":""} -b ${extend} |\\
		cut -f 2,3 > TMP/jeter2.bed

	# make title
	echo "chrom\tdelstart\tdelend\ttitle\tstart\tend" > TMP/extend.bed

	# paste to get chrom/original-start/original-end/x- start/x-end
	paste TMP/jeter1.bed TMP/jeter2.bed |\\
	    sort -T TMP -t '\t' -k1,1V -k2,2n -k3,3n --unique  >> TMP/extend.bed

	test -s TMP/extend.bed

	mv TMP/extend.bed ./

	touch versions.yml
	"""
	}


process DRAW_COVERAGE {
	tag "${row.chrom}:${row.delstart}-${row.delend} ${row.title} ${bam.name} len=${1+(row.end as int)-(row.start as int)}"
	label "process_single"
	conda "${moduleDir}/../../conda/bioinfo.01.yml"
    afterScript "rm -rf TMP"
	input:
		tuple val(meta1),path(fasta)
		tuple val(meta2),path(fai)
		tuple val(meta3),path(dict)
		tuple val(meta4),path(gtf),path(gtfidx)
		tuple val(row),val(meta),path(bam),path(bai)
	output:
		tuple val(row),path("*.pdf"),emit:pdf
		path("versions.yml"),emit:versions
	script:
		def LARGE_SV = 200_000
		def svlen = 1 + (row.end as int ) - (row.start as int)
		def gene_type = (svlen < LARGE_SV?"exon":"gene")
		def median = svlen <= LARGE_SV
	"""
	hostname 1>&2
	mkdir -p TMP

	if ${!gtf} ; then

		touch TMP/exons.R

	else

		tabix "${gtf}" "${row.chrom}:${row.start}-${row.end}" |\\
			awk -F '\t' '(\$3=="${gene_type}") {printf("%s\t%d\t%s\\n",\$1,int(\$4)-1,\$5);}' |\\
			LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\\
			bedtools merge |\\
			awk -F '\t' '{printf("rect(%d, 0,%d, -1, col=\\\"green\\\")\\n",\$2,\$3);}' > TMP/exons.R

	fi

	# extract depth with samtools
	samtools depth --threads ${task.cpus}  -a -r "${row.chrom}:${row.start}-${row.end}" "${bam}" | cut -f 2,3 > TMP/depth.txt

	# debug max and mean depth
	awk 'BEGIN{M=0;T=0.0;} {V=int(\$2);T+=V;if(V>M) M=V;} END{printf("MAX DEPTH=%d MEAN=%f\\n",M,(T/NR));}' TMP/depth.txt  1>&2


cat ${moduleDir}/plot.coverage.R |\\
	m4 -P \
	-D_LARGE_SV_=${LARGE_SV} \
	-D_CONTIG_=${row.chrom} \
	-D_CHROM_=${row.chrom} \
	-D_SAMPLE_=${meta.id} \
	-D_BAM_=${bam} \
	-D_START_=${row.start} \
	-D_END_=${row.end} \
	-D_DELSTART_=${row.delstart} \
	-D_DELEND_=${row.delend} \
	-D_MEDIAN_=${median?"TRUE":"FALSE"} \
	-D_TITLE_="${row.title}" > TMP/jeter.R

##cat TMP/jeter.R 1>&2

R --vanilla < TMP/jeter.R

mv TMP/jeter.pdf "${row.chrom}_${row.start}_${row.end}.${meta.id}.pdf"

touch versions.yml
"""
}


process MERGE_PDFS {
	tag "${meta.id}"
	label "process_single"
	conda "${moduleDir}/../../conda/ghostscript.yml"
	input:
		tuple val(meta),path("PDF/*")
	output:
		tuple val(meta),path("*.pdf"),emit:pdf
		path("versions.yml"),emit:versions
	script:
		def title = meta.id
	"""
	hostname 1>&2

	find PDF/ -type l -name "*.pdf" | LC_ALL=C sort -V > jeter.txt

	test -s jeter.txt

	gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=${title}.pdf `cat jeter.txt`

	rm jeter.txt
	touch versions.yml
	"""
	}

process ZIP_ALL {
label "process_single"
input:
	tuple val(meta),path("PDF/*")
output:
	tuple val(meta),path("all.zip"),emit:zip
	path("versions.yml"),emit:versions
script:
	def dir = (params.prefix?:"")+"archive"
	"""
hostname 1>&2
mkdir -p "${dir}"
cp -v PDF/*.pdf ${dir}/


cat << EOF > ${dir}/index.html
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
<meta http-equiv="author" content="Pierre Lindenbaum Phd ">
<title>${params.prefix?:""} Coverage</title>
<script>
var files=[
EOF

find ${dir} -type f -name "*.pdf" -printf "\\"%f\\"\\n" | paste -sd ','  >> ${dir}/index.html 

cat << EOF >> ${dir}/index.html
];
var page =0;

function goTo(dx) {
    if(files.length==0) return;
    page = page+dx;
    if(page<0) page = files.length-1;
    if(page>=files.length) page=0;
    document.getElementById("id1").src = files[page];
    document.getElementById("h1").textContent = files[page]+ " ("+(page+1)+"/"+files.length+")";
    }

 

window.addEventListener('load', (event) => {
  var frame = document.getElementById("id1");
  frame.style.height=(frame.contentWindow.document.body.scrollHeight+20)+'px';
  goTo(0);
});

</script>
</head>
<body>
<div>
    <button onclick="goTo(-1);">PREV</button>
    <span id="h1"></span>
    <button onclick="goTo(1);">NEXT</button>
</div>

<iframe id="id1" style="height:800px;width:100%;" src="blank_">
</iframe>

<div>
    <button onclick="goTo(-1);">PREV</button>
    <span>navigation</span>
    <button onclick="goTo(1);">NEXT</button>
</div>


</body>
</html>
EOF



zip -9r  "all.zip" ${dir}
rm  -f ${dir}/*.pdf
rm  -f ${dir}/*.html
touch versions.yml
"""
}

