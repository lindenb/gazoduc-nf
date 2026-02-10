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


include { dumpParams;runOnComplete         } from '../../modules/utils/functions.nf'
include { DOWNLOAD_CYTOBAND                } from '../../modules/ucsc/download.cytobands'
include { DOWNLOAD_REFGENE                 } from '../../modules/ucsc/download.refgene'
include { PREPARE_ONE_REFERENCE            } from '../../subworkflows/samtools/prepare.one.ref'

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
	 def metadata = [
		id : file(params.fasta).simpleName,
		name: file(params.fasta).simpleName
	 	]
	PREPARE_ONE_REFERENCE(
			metadata.plus(
				skip_scatter : true
				),
			Channel.fromPath(params.fasta).map{f->[[id:f.baseName],f]}
			)
  	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)

	 
	bams_ch = Channel.fromPath(params.samplesheet).
		splitCsv(sep:',',header:true).
	
	input_ch = Channel.fromPath(params.input).
		splitCsv(sep:',',header:true).
	
	ch1 = bams_ch.combine(input_ch)
		.filter{meta1,meta2->meta1.sample==meta2.sample}
		.map{meta1,meta2->meta1.plus(meta2)}
		.map{meta->{
			def contig = meta.contig;
			def start = meta.start;
			def end = meta.end;
			return [
				[id:"${contig}_${start}_${end}",contig:"${contig}",start:(start as int),end:(end as int)],
				[file(meta.bam),file(meta,bai)]
				]
			}
		.groupTuple()
		.map{meta,files->[meta,files.flatten().sort()]}

	/*
      FIND_COVERAGE_AT_LOC(
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
			PREPARE_ONE_REFERENCE.out.dict,
			bams_ch.map{meta->[file(meta.bam),file(meta.bai)]}
				.flatMap()
				.collect()
				.map{files->[[id:"bams"],files.sort()]},
			ch1.map{meta,f->meta}.unique()
			)
		ch2 = MERGE_COVERAGE_AT_LOC(file(params.sample2collection), FIND_COVERAGE_AT_LOC.out.tsv)
		*/

		
		DOWNLOAD_CYTOBAND(
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
			PREPARE_ONE_REFERENCE.out.dict
			)
		DOWNLOAD_REFGENE(
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
			PREPARE_ONE_REFERENCE.out.dict
			)

		APPLY_IGVREPORT(
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
			PREPARE_ONE_REFERENCE.out.dict,
			DOWNLOAD_CYTOBAND.out.bed,
			DOWNLOAD_REFGENE.out.bed,
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
label "process_single"
tag "${meta.id}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta4),path("BAMS/*")
	val(meta)
output:
	tuple val(meta),path("*.tsv"),emit:tsv
script:
	def prefix = task.ext.prefix?:meta.id;
"""
mkdir -p TMP
find BAMS -name "*.bam" -o -name "*.cram" | sort -V > TMP/jeter.list


jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP findallcoverageatposition -p "${meta.contig}:${meta.start}" -R ${fasta}  < TMP/jeter.list > TMP/jeter.tsv
mv TMP/jeter.tsv "${prefix}.tsv"
"""
}

process MERGE_COVERAGE_AT_LOC {
label "process_single"
tag "${key[0]}:${key[1]}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	path(sample2collection)
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
export LC_ALL=C

find COV -type l -name "*.tsv" -exec cat '{}' ';' 2> /dev/null |\\
	grep "#" -m1 > TMP/merged.tsv

find COV -type l -name "*.tsv" -exec cat '{}' ';' |\\
	grep -v "#" |\\
	sort -T TMP |\\
	uniq |\\
	sort -t \$'\\t' -k5,5  -T TMP >> TMP/merged.tsv

head TMP/merged.tsv 1>&2

echo -e "SAMPLE\tCollection" > TMP/groups.txt

if ${sample2collection.name.contains(".")}
then
	sort -T TMP -t \$'\\t' -k1,1 '${sample2collection}' --unique >> TMP/groups.txt
fi


join --header \\
	-o '1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,2.2' \\
	-t \$'\\t' \\
	-a 1 -1 5 -2 1 \\
	-e undefined  TMP/merged.tsv TMP/groups.txt > TMP/joined.tsv

mv TMP/joined.tsv TMP/merged.tsv

NROWS=`tail -n +2 TMP/merged.tsv | wc -l`

# order on the number of base that is NOT the reference
# Z is the sum of the other bases
# N0 is the number of samples with Z==0
awk -F '\t' 'BEGIN {N0=0;} (NR==1) {print;next;} {Z=0;A=int(\$16);C=int(\$17);G=int(\$18);T=int(\$19);Z=A+C+G+T; if(\$4=="A") Z-=A; if(\$4=="C") Z-=C; if(\$4=="G") Z-=G; if(\$4=="T") Z-=T; if(Z==0) {N0++;if(${max_score_0}>-1 && N0 > ${max_score_0}) next;} printf("%d\t%s\\n",Z,\$0);}'  TMP/merged.tsv > TMP/jeter1.txt
# save header 
head -n1 TMP/jeter1.txt  >  TMP/jeter2.txt

# sort on order in the first field
tail -n +2  TMP/jeter1.txt |\\
	LC_ALL=C sort -T TMP -t '\t' -k1,1nr -k6,6V | cut -f 2- >> TMP/jeter2.txt

awk -F '\t' -vT=\${NROWS} '(NR==1) {printf("%s\tPAGE\tNPAGES\\n",\$0);next;} {printf("%s\t%d\t%d\\n",\$0,1+((NR-1)/${per_page}),1+(int(T)/${per_page}));}' TMP/jeter2.txt > TMP/merged2.tsv

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
    printf("<th>Collection</th>");
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
    printf("<td><a href=\\"page%s.html\\">%s</a></td>",\$24,\$5);
    printf("<td>%s</td>",\$23);
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

tail -n +2 TMP/merged.tsv | cut -f1,24,25 > pages.tsv

awk -f TMP/jeter.awk TMP/merged.tsv > index.html
"""
}


process APPLY_IGVREPORT {
tag "${contig}:${pos} page ${page}/${page_max} N=${bams.size()}"
label "process_single"
conda "${moduleDir}/../../conda/igv-reports.yml"
//afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta4),path(cytoband)
	tuple val(meta5),path(refgene)
	tuple val(contig),val(pos),val(page),val(page_max),val(bams)
output:
	tuple val(contig),val(pos),path("page${page}.html"),emit:output
script:
	def title = contig+"_"+pos
	def refgene2 = refgene.find{it.name.endsWith(".txt.gz")}.join(" ")
	def pagei = (page as int)
	def page_maxi = (int)(page_max as double)
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
