/*

Copyright (c) 2022 Pierre Lindenbaum

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

include {moduleLoad;isHg19;isHg38;getVersionCmd} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {SIMPLE_ZIP_01} from '../../modules/utils/zip.simple.01.nf'


gazoduc = gazoduc.Gazoduc.getInstance(params)

gazoduc.build("indexcov_treshold",0.15).desc("DEL if x < 0.5+t. DUP if x > 1.5-t").menu("circos").setDouble().put()
gazoduc.build("indexcov_min_sv_length",100_000).desc("minimum SV size. Ignore events having length<x.").menu("circos").setLong().put()
gazoduc.build("indexcov_merge_length",1).desc("merge SV with distance < 'x'.").menu("circos").setLong().put()
gazoduc.build("indexcov_nsamples_cleanup",10).desc("if there is more than 'x' samples, ignore data if all samples are DEL/DUP. Ignore if < 0.").menu("circos").setInt().put()


workflow INDEXCOV_TO_CIRCOS_01 {
	take:
		meta
		reference
		bed
	main:
		version_ch  = Channel.empty()
		
		sn_ch = EXTRACT_SAMPLES(meta, bed)
		version_ch = version_ch.mix(sn_ch.version)

		one_ch = ONE_SAMPLE(meta, reference, sn_ch.filter, sn_ch.output.splitCsv(header:true,sep:'\t'))
		version_ch = version_ch.mix(one_ch.version)

		plot_ch = PLOT_CIRCOS(meta, reference, one_ch.output.collect())		
		version_ch = version_ch.mix(plot_ch.version)
		
		version_ch = MERGE_VERSION(meta, "IndexCovCircos", "IndexCov Circos",version_ch.collect())
		
	emit:
		version= version_ch
		png = plot_ch.png
		svg = plot_ch.svg
	}



process EXTRACT_SAMPLES {
executor "local"
tag "${bed.name}"
input:
	val(meta)
	path(bed)
output:
	path("samples.txt"),emit:output
	path("filter.awk"),emit:filter
	path("version.xml"),emit:version
script:
	def indexcov_treshold = meta.indexcov_treshold
	def n_cleanup = meta.indexcov_nsamples_cleanup
"""
hostname 1>&2

# number of samples
NS=`${bed.name.endsWith(".gz")?"gunzip -c ":"cat"} "${bed}" | head -n1| tr "\t" "\\n" | tail -n +4 | wc -l`

${bed.name.endsWith(".gz")?"gunzip -c ":"cat"} "${bed}" |\
	head -n 1 |\
	tr "\\t" "\\n" |\
	awk -F '\t' -vNS=\${NS} 'BEGIN{printf("column\tsample\tbed\tnsamples\\n");} (NR>3){printf("%d\t%s\t${bed.toRealPath()}\t%d\\n",NR,\$0,NS);}'  > samples.txt
test -s samples.txt


cat << 'EOF' > filter.awk
BEGIN {
	FS="\t";
	OFS="\t";
	}
/^#/	{
	next;
	}
(\$1 ~ /^(chr)?[0-9XY]+\$/)	{
	NSAMPLES = (NF - 3);
	if(${n_cleanup} > 0 && NSAMPLES >= ${n_cleanup}) {
		N_HOMDEL =0;
		N_HOMVAR =0;
		N_DUP =0;
		N_DEL =0;
		for(i=4;i<=NF;i++) {
			if( \$i < ${indexcov_treshold} ) N_HOMDEL++;
			if( \$i > (2.0 - ${indexcov_treshold} )) N_HOMVAR++;
			if( \$i > (1.5 - ${indexcov_treshold} ) && \$i < (1.5 + ${indexcov_treshold}  ) ) N_DUP++;
			if( \$i < (0.5 + ${indexcov_treshold} )  && \$i > (0.5 - ${indexcov_treshold})) N_DEL++;
			}
		if(N_HOMDEL >= NSAMPLES) next;
		if(N_HOMVAR >= NSAMPLES) next;
		if(N_DUP >= NSAMPLES) next;
		if(N_DEL >= NSAMPLES) next;
		}
	print;
	}
EOF

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">Extract samples and their index in indexcov file</entry>
        <entry key="bed">${bed}</entry>
        <entry key="versions">${getVersionCmd("awk")}</entry>
</properties>
EOF
"""
}


process ONE_SAMPLE {
tag "${row.sample}"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
	path(filter)
	val(row)
output:
	path("data.conf"),emit:output
	path("version.xml"),emit:version
script:
	def indexcov_treshold = meta.indexcov_treshold?:0.2
	def nsamples = (row.nsamples as double)
	def R0=0.1
	def R1=0.9
	def DR=(R1-R0)/nsamples
	def r0= R0 + (DR*((row.column as int)-1))
	def r1 = r0 + DR
	def min_size = meta.indexcov_min_sv_length
	def merge_d = meta.indexcov_merge_length
	def url = isHg38(reference)?"https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=[chr]%3A[start]%2D[end]":""
"""
hostname 1>&2
set -o pipefail
${moduleLoad("bedtools")}
mkdir -p TMP

${row.bed.endsWith(".gz")?"gunzip -c ":"cat"} "${row.bed}" |\
	awk -F '\t' -f "${filter}" |\
	cut -d '\t' -f '1,2,3,${row.column}' > TMP/tmp1.bed

awk -F '\t' '( \$4 <= (0.5 + ${indexcov_treshold} ) )'  TMP/tmp1.bed  |\
	sort -T TMP -t '\t' -k1,1 -k2,2n  |\
	bedtools merge -d '${merge_d}' -o mean -c 4|\
	awk '{L=int(\$3)-int(\$2)+1;if(L< ${min_size}) next; C=(\$4 <= ${indexcov_treshold} ?"red":"orange");printf("%s\t%s\t%s\tfill_color=%s,stroke_color=%s\\n",\$1,\$2,\$3,C,C);}' > TMP/del.bed

awk -F '\t' '( \$4 >= (1.5 - ${indexcov_treshold} ) )'  TMP/tmp1.bed  |\
	sort -T TMP -t '\t' -k1,1 -k2,2n  |\
	bedtools merge -d '${merge_d}'  -o mean -c 4 |\
	awk '{L=int(\$3)-int(\$2)+1;if(L< ${min_size}) next;C=(\$4 >= 2.0 - ${indexcov_treshold} ?"blue":"cyan"); printf("%s\t%s\t%s\tfill_color=%s,stroke_color=%s\\n",\$1,\$2,\$3,C,C);}' > TMP/dup.bed
	

cat TMP/del.bed TMP/dup.bed | sed 's/^chr/hs/' | sort -T TMP -t '\t' -k1,1 -k2,2n > "${row.sample}.circos"


cat << EOF > data.conf
<highlight>
file= \${PWD}/${row.sample}.circos
r0=${r0}r
r1=${r1}r
${url.isEmpty()?"":"url=${url}"}
stroke_thickness = 2
</highlight>
EOF



###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">build sample circors conf</entry>
        <entry key="bed">${row.bed}</entry>
        <entry key="sample">${row.sample}</entry>
        <entry key="versions">${getVersionCmd("bedtools awk")}</entry>
</properties>
EOF
"""
}


process PLOT_CIRCOS {
tag "N=${L.size()}"
input:
	val(meta)
	val(reference)
	val(L)
output:
	path("${meta.prefix?:""}circos.svg"),emit:svg
	path("${meta.prefix?:""}circos.png"),emit:png
	path("version.xml"),emit:version
script:
	def karyotype = isHg38(reference)?"data/karyotype/karyotype.human.hg38.txt":(isHg19(reference)?"data/karyotype/karyotype.human.hg19.txt":"")
"""
hostname 1>&2
${moduleLoad("circos")}


cat << EOF > jeter.conf
karyotype = ${karyotype}
<ideogram>
<spacing>
default = 0.005r
</spacing>
radius = 0.9r
thickness = 20p
fill      = yes
</ideogram>
<image>
<<include etc/image.conf>>
</image>

<highlights>
EOF

cat ${L.join(" ")} >> jeter.conf

cat << EOF >> jeter.conf
</highlights>
<<include etc/colors_fonts_patterns.conf>>
<colors>
reda = 200,0,0,0.3
bluea=0,0,200,0.3
</colors>
<<include etc/housekeeping.conf>>
EOF



circos --conf jeter.conf  -outputdir . 

mv circos.svg "${meta.prefix?:""}circos.svg"
mv circos.png "${meta.prefix?:""}circos.png"


###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">build sample circors conf</entry>
        <entry key="versions">${getVersionCmd("awk")}</entry>
</properties>
EOF
"""
}

