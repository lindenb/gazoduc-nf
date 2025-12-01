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
include {GHOSTSCRIPT_MERGE                } from '../../modules/gs/merge'
include {PDF_NAVIGATION                   } from '../../modules/pdf/navigation'

workflow PLOT_COVERAGE_01 {
	take:
		metadata
		fasta
		fai
		dict
		gtf
		bed // [meta,bed]
		bams //[meta,bam,bai]
	main:
		versions = Channel.empty()
		EXTEND_BED(fasta,fai,dict, bed)
		versions = versions.mix(EXTEND_BED.out.versions)

		c1_ch = EXTEND_BED.out.bed.
			map{_meta,f->f}.
			splitCsv(header: true,sep:'\t',strip:true)

		DRAW_COVERAGE(fasta,fai,dict,gtf,c1_ch.combine(bams))
		versions = versions.mix(DRAW_COVERAGE.out.versions)

		GHOSTSCRIPT_MERGE(
			DRAW_COVERAGE.out.pdf
				.map{T->[[id:T[0].chrom+"_"+T[0].delstart+"_"+T[0].delend],T[1]]}
				.groupTuple()
			)
		versions = versions.mix(GHOSTSCRIPT_MERGE.out.versions)


		PDF_NAVIGATION(GHOSTSCRIPT_MERGE.out.pdf.map{_meta,pdf->pdf}.collect().map{[metadata,it]})
		versions = versions.mix(PDF_NAVIGATION.out.versions)
	emit:
		zip =PDF_NAVIGATION.out.zip
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

	{
	L=int(\$3)-int(\$2);
	if(M>0 && L>M) next;
	if(N>0 && L<N) next;
	printf("%s\t%s\t%s\t%s\\n",\$1,\$2,\$3, (NF==3 || \$4==""?".":\$4) );
	}
EOF

	grep -vE '^(chrom|browser|track|#)'  "${bed}" |\\
	jvarkit  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP bedrenamechr -R "${fasta}" -c 1 | \\
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
	stub:
	"""
	touch versions.yml
	echo "chrom\tdelstart\tdelend\ttitle\tstart\tend" >  extend.bed
	awk -F '\t' '{printf("%s\t%s\t%s\ttitle\t%s\t%s\\n",\$1,\$2,\$3,\$2,\$3);}' ${bed} >> extend.bed
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
stub:
"""
touch versions.yml "${row.chrom}_${row.start}_${row.end}.${meta.id}.pdf"
"""
}

