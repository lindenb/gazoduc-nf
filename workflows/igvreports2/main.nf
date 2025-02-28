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
nextflow.enable.dsl=2


include {dumpParams;runOnComplete} from '../../modules/utils/functions.nf'
//include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {DOWNLOAD_CYTOBAND} from '../../modules/ucsc/download.cytobands'
include {DOWNLOAD_REFGENE} from '../../modules/ucsc/download.refgene'

if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}

def toLoc(def row) {
	if(row.size()==1) {
		}
	else if(row.size()==2) {
		return row;
		}
	}

workflow {
	 genome = Channel.of(file(params.fasta), file(params.fai), file(params.dict)).collect()
	 bams_ch = Channel.fromPath(params.bams).
		splitText().
		map{it.trim()}.
		collate(50)		

	Channel.fromPath(params.positions).
		splitCsv(sep:'\t',header:false).
		map{T->{
			switch(T.size()) {
				case 2: T[1] = [T[0],(T[1] as int)];
				case 1: 
					int colon =T[0].lastIndexOf(":");
					if(colon==-1) throw new IllegalArgumentException("cannot find ':' in "+T[0]);
					int pos = (T[0].substring(colon+1) as int)
					return [T[0].substring(0,colon), pos];
				default: throw new IllegalArgumentException("splitCsv from position");
				}
			}}.combine(bams_ch).
			map{[it[0],it[1],it[2 ..< it.size()].sort()]}.
			set{contig_pos_bams_ch}


         	covpos_ch = FIND_COVERAGE_AT_LOC(genome,contig_pos_bams_ch)
		ch1 = covpos_ch.map{[ [it[0],it[1]], it[2]]}.groupTuple()
		ch2 = MERGE_COVERAGE_AT_LOC(ch1)
		
		ch4 = ch2.output.
			map{[it[1],it[0]]}.
			splitCsv(sep:'\t',header:false).
			map{[ [it[1][0] /* contig */ ,it[1][1] /* pos */,it[0][1] /* page */,it[0][2] /*pages*/], it[0][0] /* bam */] }.
			groupTuple().
			map{[it[0][0],it[0][1],it[0][2],it[0][3],it[1]]}

		cyto_ch = DOWNLOAD_CYTOBAND(genome)

		refgene_ch = DOWNLOAD_REFGENE(genome)

		report_ch = APPPLY_IGVREPORT(
			genome,
			cyto_ch.output,
			refgene_ch.output,
			ch4
			)
		
		ch5 = report_ch.output.map{[[it[0],it[1]],it[2]]}.
			mix(ch2.output.map{[[it[0][0],it[0][1]],it[2]]}).
			groupTuple().
			map{[it[0][0],it[0][1],it[1].sort()]}

		ch6 = MAKE_DIRECTORY(ch5)

		ZIPIT(ch6.output.collect(),ch6.tsv.collect())

	}

runOnComplete(workflow)


process FIND_COVERAGE_AT_LOC {
label "process_quick"
tag "${contig}:${pos} N=${bams.size()}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	path(genome)
	tuple val(contig),val(pos),val(bams)
output:
	tuple val(contig),val(pos),path("*.tsv"),emit:output
script:
	def fasta = genome.find{it.name.endsWith("a")}
"""
set -o pipefail
mkdir -p TMP
cat << '__EOF__' > TMP/jeter.list
${bams.join("\n")}
__EOF__

MD5=`cat TMP/jeter.list | md5sum | cut -d ' ' -f1`

jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP findallcoverageatposition -p "${contig}:${pos}" -R ${fasta}  < TMP/jeter.list > TMP/jeter.tsv
mv TMP/jeter.tsv "coverage.\${MD5}.tsv"
"""
}

process MERGE_COVERAGE_AT_LOC {
label "process_quick"
tag "${key[0]}:${key[1]}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	tuple val(key),path("COV/*")
output:
	tuple val(key),path("pages.tsv"),path("index.html"),emit:output
script:
	def contig = key[0];
	def pos = key[1];
	def per_page=20
	def max_score_0 = params.max_score_0?:-1;
"""
mkdir -p TMP
find COV -type l -name "*.tsv" -exec cat '{}' ';' |\\
	LC_ALL=C sort -T TMP |\\
	uniq > TMP/merged.tsv

head TMP/merged.tsv 1>&2

NROWS=`tail -n +2 TMP/merged.tsv | wc -l`

# order on the number of base that is NOT the reference
# Z is the sum of the other bases
# N0 is the number of samples with Z==0
awk -F '\t' 'BEGIN {N0=0;} (NR==1) {print;next;} {Z=0;A=int(\$16);C=int(\$17);G=int(\$18);T=int(\$19);Z=A+C+G+T; if(\$4=="A") Z-=A; if(\$4=="C") Z-=C; if(\$4=="G") Z-=G; if(\$4=="T") Z-=T; if(Z==0) {N0++;if(${max_score_0}>=0 ||N0>{max_score_0}) next;} printf("%d\t%s\\n",Z,\$0);}'  TMP/merged.tsv > TMP/jeter1.txt
# save header 
head -n1 TMP/jeter1.txt  >  TMP/jeter2.txt

# sort on order in the first field
tail -n +2  TMP/jeter1.txt |\\
	LC_ALL=C sort -T TMP -t '\t' -k1,1nr -k6,6V | cut -f 2- >> TMP/jeter2.txt

awk -F '\t' -vT=\${NROWS} '(NR==1) {printf("%s\tPAGE\tNPAGES\\n",\$0);next;} {printf("%s\t%d\t%s\\n",\$0,1+((NR-1)/${per_page}),1+(int(T)/${per_page}));}' TMP/jeter2.txt > TMP/merged2.tsv

mv TMP/merged2.tsv TMP/merged.tsv


cat << '__EOF__' > TMP/jeter.awk
BEGIN {
	printf("<!DOCTYPE html>\\n<html><head><meta charset=\\"UTF-8\\"><title>${contig}:${pos}</title></head><body><h1>${contig}:${pos}</h1><table border=\\"1\\">");
	}

NR == 1 {
    printf("<thead><tr>");
    printf("<th>Position</th>");
    printf("<th>REF</th>");
    printf("<th>Sample</th>");
    printf("<th>DEPTH</th>");
    printf("<th>A</th>");
    printf("<th>C</th>");
    printf("<th>G</th>");
    printf("<th>T</th>");
    printf("</tr></thead><tbody>\\n");
}

NR > 1 {
    printf("<tr>");
    printf("<td>%s:%s</td>",\$2,\$3);
    printf("<td>%s</td>",\$4);
    printf("<td><a href=\\"page%s.html\\">%s</a></td>",\$23,\$5);
    printf("<td>%s</td>",\$6);
    printf("<td>%s</td>",\$16);
    printf("<td>%s</td>",\$17);
    printf("<td>%s</td>",\$18);
    printf("<td>%s</td>",\$19);
    printf("<tr>\\n");
}

END {
	printf("</tbody></table></body></html>\\n");
	}
__EOF__

tail -n +2 TMP/merged.tsv | cut -f1,23,24 > pages.tsv

awk -f TMP/jeter.awk TMP/merged.tsv > index.html
"""
}


process APPPLY_IGVREPORT {
tag "${contig}:${pos} page ${page}/${page_max} N=${bams.size()}"
label "process_quick"
conda "${moduleDir}/../../conda/igv-reports.yml"
//afterScript "rm -rf TMP"
input:
	path(genome)
	path(cytoband)
	path(refgene)
	tuple val(contig),val(pos),val(page),val(page_max),val(bams)
output:
	tuple val(contig),val(pos),path("page${page}.html"),emit:output
script:
	def title = contig+"_"+pos
	def fasta = genome.find{it.name.endsWith("a")}
	def refgene2 = refgene.find{it.name.endsWith(".txt.gz")}.join(" ")
	def pagei = (page as int)
	def page_maxi = (page_max as int)
        def navigation= "<a href=\"page"+(pagei==1?page_maxi:pagei-1)+".html\">[Previous]</a> " +
		"<a href=\"index.html\">[Index]</a> " +
		"<a href=\"page"+(pagei+1>page_maxi?1:pagei+1)+".html\">[Next]</a> "
"""
hostname 1>&2
mkdir -p TMP

echo "${contig}\t${(pos as int)-1}\t${pos}" > TMP/jeter.bed

create_report TMP/jeter.bed  "${fasta}" \\
	--ideogram "${cytoband}" \\
	--tracks ${bams.join(" ")} ${refgene2} \\
	--output TMP/jeter.html

cat << '__EOF__' > TMP/jeter.xsl
<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output method="html" encoding="UTF-8"/>
  
  <xsl:template match="/">
    <xsl:apply-templates/>
  </xsl:template>

  <xsl:template match="title[name(..)='head']">
  <title>${contig}:${pos}</title>
  </xsl:template>
   
  <xsl:template match="body[name(..)='html']">
  <body>
  <h1>${contig}:${pos}</h1>
  <div>${navigation}</div>
  <xsl:apply-templates/>
  <div>${navigation}</div>
  </body>
  </xsl:template>

  <xsl:template match="@*|node()">
        <xsl:copy>
          <xsl:copy-of select="@*"/>
          <xsl:apply-templates/>
        </xsl:copy>
  </xsl:template>

</xsl:stylesheet>
__EOF__

xsltproc --html  TMP/jeter.xsl TMP/jeter.html  > TMP/jeter2.html

mv -v TMP/jeter2.html ./page${page}.html
"""
}
process MAKE_DIRECTORY {
tag "${contig}_${pos}"
executor "local"
input:
	tuple val(contig),val(pos),path("HTML/*")
output:
	path("${contig}_${pos}"),emit:output
	path("${contig}_${pos}.tsv"),emit:tsv
script:
	
"""
mkdir -p "${contig}_${pos}"
cp  HTML/*.html "${contig}_${pos}/"
echo "${contig}\t${pos}\t<tr><td><a href='${contig}_${pos}/index.html'>${contig}:${pos}</a></td></tr>" > ${contig}_${pos}.tsv
"""
}



process ZIPIT {
executor "local"
input:
	path("DIR/*")
	path("HTML/*")
output:
	path("archive.zip"),emit:output
script:
"""
mkdir -p archive

cat << EOF >> archive/index.html
<!DOCTYPE html>
<html><head><meta charset=\\"UTF-8\\"><title>${params.prefix?:""}IGV reports</title></head>
<body>
<table>
 <tbody>
EOF

cp -vr DIR/* archive/

find HTML -name "*.tsv" -exec cat '{}' ';' |\\
	sort -t '\t' -k1,1V -k2,2n | cut -f3- >> archive/index.html

echo "</tbody></table></body></html>" >> archive/index.html

zip -9r archive.zip archive
"""
}
