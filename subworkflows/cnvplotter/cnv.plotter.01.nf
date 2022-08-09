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
include {SAMTOOLS_SAMPLES01} from '../../modules/samtools/samtools.samples.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {moduleLoad;isBlank;isHg38;isHg19} from '../../modules/utils/functions.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'

workflow CNV_PLOTTER_01 {
	take:
		meta	
		reference
		vcf
		bams
	main:
		version_ch = Channel.empty()

		ch1_ch = SAMTOOLS_SAMPLES01([:],reference,bams)
		version_ch = version_ch.mix(ch1_ch.version)

		compile_ch = COMPILE_VCF_PARSER(meta)
		version_ch = version_ch.mix(compile_ch.version)
		
		splitctx_ch = SPLIT_VARIANTS(meta,vcf,compile_ch.jar)		
		version_ch = version_ch.mix(splitctx_ch.version)

		ch2_ch = splitctx_ch.output.splitCsv(header:true,sep:'\t').
			combine(ch1_ch.output).
			map{T->T[0].plus("bams":T[1])}


		plot_ch = PLOT_CNV(meta, reference, ch2_ch)
		version_ch = version_ch.mix(plot_ch.version)

		merge_ch = MERGE_PLOTS(meta,plot_ch.output.collect())		
		version_ch = version_ch.mix(merge_ch.version)

		version_ch = MERGE_VERSION(meta, "CNVPlot", "CNV Plotter", version_ch.collect())
	emit:
		version = version_ch
		zip = merge_ch.zip
	}


process COMPILE_VCF_PARSER {
executor "local"
afterScript "rm -rf TMP"
input:
	val(meta)
output:
	path("minikit.jar"),emit:jar
	path("version.xml"),emit:version

script:
"""
${moduleLoad("jvarkit")}

mkdir TMP

cat << "__EOF__" > TMP/Minikit.java
import java.io.BufferedReader;
import java.util.regex.*;
import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.stream.*;
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.*;
import htsjdk.variant.vcf.*;



public class Minikit {

private final int minLenOnReference = ${meta.minCnvLength?:"1_000"};
private final int maxLenOnReference = ${meta.maxCnvLength?:"10_000_000"};
private final String prefix="${meta.prefix?:""}";
private final int max_controls =  ${meta.max_controls?:"50"};

private boolean hasCNV(Genotype g) {	
	return g.isHet() || g.isHomVar();
	}

private void instanceMain(final String args[]) {
	try {
		String vcf="";
		int optind=0;
		while(optind < args.length) {
			if(args[optind].equals("--vcf") && optind+1< args.length) {
				optind++;
				vcf = args[optind];
				}
			else if(args[optind].equals("--")) {
				optind++;
				break;
				}
			else if(args[optind].startsWith("-")) {
				System.err.println("unknown option "+args[optind]);
				System.exit(-1);
				}
			optind++;
			}
		if(optind != args.length) {
			System.err.println("illegal number of arguments");
			System.exit(-1);
			}
		PrintStream out = System.out;
		out.print("prefix");
		out.print("\t");
		out.print("interval");
		out.print("\t");
		out.print("vcf");
		out.print("\t");
		out.print("cases");
		out.print("\t");
		out.print("controls");
		out.print("\t");
		out.print("html");
		out.println();
		try(VCFIterator r = new VCFIteratorBuilder().open(System.in)) {
			long variant_id =0L;
			final VCFHeader h = r.getHeader();
			while(r.hasNext()) {
				final VariantContext ctx = r.next();
				int len = ctx.getLengthOnReference();
				if(len < minLenOnReference || len > maxLenOnReference) continue;
				final String svType = ctx.getAttributeAsString("SVTYPE","");
				if(svType==null || svType.isEmpty() || svType.equals("INV")) continue;
				final Set<String> affected = ctx.getGenotypes().stream().
					filter(G->hasCNV(G)).
					map(G->G.getSampleName()).
					collect(Collectors.toSet());
				if(affected.isEmpty()) continue;
				List<String> unaffected = ctx.getGenotypes().stream().
					filter(G->!hasCNV(G)).
					map(G->G.getSampleName()).
					collect(Collectors.toList());
				final StringBuilder sb = new StringBuilder();
				sb.append("Variant location is <b>"+ctx.getContig()+":"+ctx.getStart()+"-"+ctx.getEnd()+". ");
				sb.append("SVType is <b>"+svType+"</b>. ");
					

				Collections.shuffle(unaffected);
				if(unaffected.size()  > max_controls)  unaffected = unaffected.subList(0,max_controls);
				variant_id++;
				out.print(this.prefix + ctx.getContig()+"_"+ctx.getStart()+"_"+ctx.getEnd()+"_"+svType+"_id"+variant_id);
				out.print("\t");
				out.print(ctx.getContig()+":"+ctx.getStart()+"-"+ctx.getEnd());
				out.print("\t");
				out.print(vcf);
				out.print("\t");
				out.print(String.join(",",affected));
				out.print("\t");
				out.print(String.join(",",unaffected));
				out.print("\t");
				out.print(sb.toString());
				out.println();
			}// while r.hasNext
		out.flush();
		}//try
	    }
	catch(final Throwable err ) {
		err.printStackTrace();
		System.exit(-1);
	}
}
public static void main(final String[] args) {
	new Minikit().instanceMain(args);
	}
}
__EOF__


cat << EOF > TMP/tmp.mf
Manifest-Version: 1.0
Main-Class: Minikit
EOF


javac -cp \${JVARKIT_DIST}/coverageplotter.jar -d TMP -sourcepath TMP TMP/Minikit.java
jar cfm minikit.jar TMP/tmp.mf -C TMP .

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">compile minikit</entry>
	<entry key="javac.version">\$(javac -version 2>&1)</entry>
</properties>
EOF
"""
stub:
"""
touch "minikit.jar"
echo "<properties/>" > version.xml
"""
}


process SPLIT_VARIANTS {
executor "local"
afterScript "rm -rf TMP"
input:
      	val(meta)
	val(vcf)
	val(minikit)
output:
       	path("variants.tsv"),emit:output
        path("version.xml"),emit:version
script:
	def extra_filter = meta.extra_vcf_filter?:""
"""
hostname 1>&2
${moduleLoad("bcftools jvarkit")}
set -o pipefail

bcftools view "${vcf}" |\
	${isBlank(extra_filter)?"":"${extra_filter} |"} \
	java -cp  \${JVARKIT_DIST}/coverageplotter.jar:${minikit} Minikit --vcf "${vcf}" > variants.tsv



###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">extract variants</entry>
	<entry key="vcf">${vcf}</entry>
</properties>
EOF
"""
}

process PLOT_CNV {
tag "${row.prefix} ${row.interval}"
memory "3g"
input:
	val(meta)
	val(reference)
	val(row)
output:
	tuple path("${row.prefix}out.html"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("jvarkit")}
set -o pipefail
mkdir TMP

sort -T TMP -t '\t' -k1,1 "${row.bams}" > TMP/samples.bams.tsv

echo "${row.cases}" | tr "," "\\n" | sort | uniq > TMP/cases.txt
echo "${row.controls}" | tr "," "\\n" | sort | uniq > TMP/controls.txt

join -t '\t' -1 1 -2 1 -o "2.2" TMP/cases.txt TMP/samples.bams.tsv > TMP/cases.bams.list
join -t '\t' -1 1 -2 1 -o "2.2" TMP/controls.txt TMP/samples.bams.tsv | shuf | head -n 10 > TMP/controls.bams.list
cat TMP/cases.bams.list TMP/controls.bams.list  | sort | uniq > TMP/all.bams.list

java  -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP  -jar ${JVARKIT_DIST}/coverageplotter.jar \
	-R "${reference}" --region "${row.interval}" `cat TMP/all.bams.list` > TMP/jeter.html

cat << EOF > TMP/jeter.xsl
<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"  version="1.0">

  <xsl:template match="div[@id='__PLACEHOLDER__']">
   <div>
   <h3>Description</h3>
   <div>${row.html}</div>
   </div>
  </xsl:template>

  <xsl:template match="@*|node()">
    <xsl:copy>
      <xsl:apply-templates select="@*|node()" />
    </xsl:copy>
  </xsl:template>

</xsl:stylesheet>
EOF


xsltproc TMP/jeter.xsl TMP/jeter.html >  "${row.prefix}out.html"

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">plot CNV</entry>
	<entry key="interval">${row.interval}</entry>
</properties>
EOF
"""
}

process MERGE_PLOTS {
tag "N=${L.size()}"
input:
	val(meta)
	val(L)
output:
	path("${meta.prefix?:""}plots.zip"),emit:zip
	path("version.xml"),emit:version
script:
	def prefix = "${meta.prefix?:""}all"
	def copy = L.collect{T->"ln -s '"+T+"' '"+prefix+"/'"}.join("\n")
	def pages = L.collect{T->"\""+file(T).name+"\""}.sort().join(",")
"""
hostname 1>&2

mkdir "${prefix}" 

${copy}

cat << EOF > ${prefix}/index.html
<html>
<head>
<title>${meta.prefix?:""} CNV plotter (N=${L.size()})</title>
<script>
var files=[
${pages}
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
</body>
</html>
EOF

zip -r -9 "${meta.prefix?:""}plots.zip" "${prefix}"

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">Zip SVG plots</entry>
</properties>
EOF
"""
}

