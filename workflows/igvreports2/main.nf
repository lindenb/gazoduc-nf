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
			}}.combine(bams_ch).set{contig_pos_bams_ch}

	

         	covpos_ch = FIND_COVERAGE_AT_LOC(genome,contig_pos_bams_ch)
		cyto_ch = DOWNLOAD_CYTOBAND(genome)

		refgene_ch = DOWNLOAD_REFGENE(genome)


		
		report_ch = IGVREPORT(genome, cyto_ch.output, refgene_ch.output, prepare_ch.output.splitCsv(header:true,sep:'\t') )
                version_ch = version_ch.mix(report_ch.version)
		
		ZIPIT(report_ch.output.collect())
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
mv TMP/jeter.tsv covarage.${MD5}.tsv
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
	tuple val(key),path("pages.tsv"),path(index.html"),emit:pages
script:
	def contig = key[0];
	def pos = key[1];
	def per_page=20
"""
mkdir -p TMP
find COV -type l -name "*.tsv" -exec cat '{}' ';' |\\
	LC_ALL=C sort -T TMP |\\
	uniq > TMP/merged.tsv

NROWS=`tail -n +2 TMP/merged.tsv | wc -l`

# order on the number of base that is NOT the reference
awk -F '\t' '(NR==1) {printf("#ORDER\t%s\\n",\$0);next;} {A=int(\$16);C=int(\$17);G=int(\$18);T=int(\$19);Z=A+C+G=T; if($4=="A") Z-=A; if(\$4=="C") Z-=C; if(\$4=="G") Z-=G; if(\$4=="T") Z-=T; printf("%s\t%s\\n",Z,$0);}'  TMP/merged.tsv |\\
	LC_ALL=C sort -T TMP -t '\t' -k1,1nr -k6,6V |\\
	cut -f2- |\\
	awk -F '\t' -vT=\${NROWS} '(NR==1) {printf("%s\tPAGE\tNPAGES\\n",\$0);next;} {printf("%s\t%d\t%s\\n",1+((NR-1)/${per_page}),1+(int(T)/${per_page});}' > TMP/merged2.tsv
mv TMP/merge2d.tsv TMP/merged.tsv


cat << __EOF__ > TMP/jeter.awk
BEGIN {
	printf("<!DOCTYPE html>\n"<html><head><meta charset=\"UTF-8\"><title>${contig}:${pos}</title></head><body><h1>${contig}:${pos}</h1><table border=\"1\"><thead>");
	}

NR == 1 {
    print "<thead><tr>";
    prinf("<th>Position</th>");
    prinf("<th>REF</th>");
    prinf("<th>Sample</th>");
    prinf("<th>DEPTH</th>");
    prinf("<th>A</th>");
    prinf("<th>C</th>");
    prinf("<th>G</th>");
    prinf("<th>T</th>");
    print "</tr></thead><tbody>\\n");
}

NR > 1 {
    printf("<tr>");
    prinf("<th>%s:%s</th>",\$2,\$3);
    prinf("<th>%s</th>",\$4);
    prinf("<th><a href=\"page%s.html\">%s</a></th>",\$23,\$5);
    prinf("<th>%s</th>",\$6);
    prinf("<th>%s</th>",\$16);
    prinf("<th>%s</th>",\$17);
    prinf("<th>%s</th>",\$18);
    prinf("<th>%s</th>",\$19);
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
tag "${row.title}"
afterScript "rm -rf TMP"
input:
	path(genome)
	path(cytoband)
	path(refgene)
	tuple val(contig),val(start),val(end),val(page),val(page_max),val(tracks)
output:
	tuple val(contig+"_"+start+(start==end?"":"_"+end)),path("*.html"),emit:output
script:
	def title = contig+"_"+start+(start==end?"":"_"+end);
	def fasta = genome.find{it.name.endsWith("a")}
	def refgene2 = refgene.find{it.name.endsWith(".txt.gz")}.join(" ")
"""
hostname 1>&2
mkdir -p TMP

create_report ${row.vcf.isEmpty()?row.bedpe:row.vcf}  ${fasta} \\
	--ideogram "${cytoband}" \\
	--tracks ${tracks} ${refgene2} \\
	--output TMP/jeter.html

cat << EOF > TMP/jeter.xsl
<xsl:template match="title">
<title>${title}</title>
</xsl:template>

<xsl:template match="body">
<body>
<div>
<a href="">Previous page</a>
<a href="">Next page</a>
</div>
<xsl:apply-templates select="*|text()">
<div>
</div>
<a href="">Previous page</a>
<a href="">Next page</a>
</body>
</xsl:template>

<xsl:template match="body">
</xsl:template>

EOF


mv -v "TMP/${row.title}.html" ./
"""
}

process ZIPIT {
input:
	path(htmls)
output:
	path("archive.zip"),emit:output
script:
"""
zip -9 -j archive.zip ${htmls}
"""
}
