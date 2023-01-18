/*

Copyright (c) 2023 Pierre Lindenbaum

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
include {DOWNLOAD_GNOMAD_SV_01} from '../../modules/gnomad/download.gnomad.sv.01.nf'
include {DOWNLOAD_DGV_01} from '../../modules/dgv/download.dgv.01.nf'
include {DOWNLOAD_GFF3_01} from '../../modules/gff3/download.gff3.01.nf'

def gazoduc = gazoduc.Gazoduc.getInstance(params);

gazoduc.make("max_cases",100_000).
	description("max number of cases per plot").
	setInt().
	put()

gazoduc.make("max_controls",10).
	description("max number of controls per plot").
	setInt().
	put()

gazoduc.make("minCnvLength",1).
        description("min CNV Length").
        setInt().
        put()

gazoduc.make("maxCnvLength",250_000_000).
        description("min CNV Length").
        setInt().
        put()



workflow CNV_PLOTTER_01 {
	take:
		meta	
		reference
		vcf
		bams
		excludeids
	main:
		version_ch = Channel.empty()

		ch1_ch = SAMTOOLS_SAMPLES01([:],reference,bams)
		version_ch = version_ch.mix(ch1_ch.version)


		merge_ch = Channel.empty()
		gnomad_ch = DOWNLOAD_GNOMAD_SV_01(meta,reference)
		version_ch = version_ch.mix(gnomad_ch.version)
		merge_ch = merge_ch.mix(gnomad_ch.bed)

		dgv_ch = DOWNLOAD_DGV_01(meta,reference)
		version_ch = version_ch.mix(dgv_ch.version)
		merge_ch = merge_ch.mix(dgv_ch.bed)

		known_ch = MERGE_KNOWN(meta,merge_ch.collect())
		version_ch = version_ch.mix(known_ch.version)

		gff3_ch = DOWNLOAD_GFF3_01(meta.plus(with_tabix:"true"),reference)
		version_ch = version_ch.mix(gff3_ch.version)

		compile_ch = COMPILE_VCF_PARSER(meta,reference)
		version_ch = version_ch.mix(compile_ch.version)
		
		splitctx_ch = SPLIT_VARIANTS(meta,vcf,excludeids, compile_ch.jar)		
		version_ch = version_ch.mix(splitctx_ch.version)

		ch2_ch = splitctx_ch.output.splitCsv(header:true,sep:'\t').
			combine(ch1_ch.output).
			map{T->T[0].plus([
				"bams":T[1],
				"max_cases":(meta.max_cases?:1000000),
				"max_controls":(meta.max_controls?:10)
				])}

		ch3_ch = ch2_ch.filter{T->!T.svtype.equals("INV")}
		ch4_ch = ch2_ch.filter{T->T.svtype.equals("INV")}
	
		plot_ch = PLOT_CNV(meta, reference, known_ch.bed, gff3_ch.gff3 , ch3_ch)
		version_ch = version_ch.mix(plot_ch.version)


		plotinv_ch = PLOT_INV(meta, reference, known_ch.bed, gff3_ch.gff3 , ch4_ch)
		version_ch = version_ch.mix(plotinv_ch.version)


		all_html = plot_ch.output.concat(plotinv_ch.output).collect()

		merge_ch = MERGE_PLOTS(meta,splitctx_ch.output,all_html)
		version_ch = version_ch.mix(merge_ch.version)

		gifch = MAKEGIF(meta,all_html)
		version_ch = version_ch.mix(gifch.version)

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
	val(reference)
output:
	path("minikit.jar"),emit:jar
	path("version.xml"),emit:version

script:
"""
${moduleLoad("jvarkit")}

mkdir TMP

cat << "__EOF__" > TMP/Minikit.java
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
import java.math.*;
import java.security.MessageDigest;
import javax.xml.*;
import javax.xml.stream.*;


public class Minikit {
private static final String RDF="http://www.w3.org/1999/02/22-rdf-syntax-ns#";
private static final String U1087="https://umr1087.univ-nantes.fr/";

private final int minLenOnReference = ${meta.minCnvLength?:"1"};
private final int maxLenOnReference = ${meta.maxCnvLength?:"250_000_000"};
private final String prefix="${meta.prefix?:""}";
private final int max_controls =  ${meta.max_controls?:"50"};

private static String  md5(String s) {
	MessageDigest md;
	 try {
		 md = MessageDigest.getInstance("MD5");
	 } catch (final Exception err) {
		throw new RuntimeException(err);
	 	}
	md.update(s.getBytes());
	return new BigInteger(1,md.digest()).toString(16);	
	}

private boolean hasCNV(Genotype g) {	
	return g.isHet() || g.isHomVar();
	}

private String getVariantId(VariantContext ctx) throws Exception {
	final String svType = ctx.getAttributeAsString("SVTYPE","undefined");
	final String ref = "${file(reference).getSimpleName()}";
	return md5(ref+":"+ctx.getContig()+":"+ctx.getStart()+"-"+ctx.getEnd()+":"+svType);
	}

private String rdf(VariantContext ctx) throws Exception {
	final String svType = ctx.getAttributeAsString("SVTYPE","undefined");
	final String ref = "${file(reference).getSimpleName()}";
	final XMLOutputFactory xof = XMLOutputFactory.newFactory();
	final StringWriter sw = new StringWriter();
	final XMLStreamWriter w = xof.createXMLStreamWriter(sw);
	final String id = getVariantId(ctx); 
	w.writeStartElement("rdf","RDF",RDF);
	w.writeNamespace("rdf", RDF);
	w.writeNamespace("u", U1087);
	w.writeAttribute(XMLConstants.XML_NS_PREFIX,XMLConstants.XML_NS_URI,"base",U1087);
	w.writeStartElement("u","Variant",U1087);
	w.writeAttribute("rdf", RDF, "about",id);


	w.writeStartElement("u", "build", U1087);
	w.writeCharacters(ref);
	w.writeEndElement();
	
	w.writeStartElement("u", "contig", U1087);
	w.writeCharacters(ctx.getContig());
	w.writeEndElement();
	
	w.writeStartElement("u", "start", U1087);
	w.writeCharacters(String.valueOf(ctx.getStart()));
	w.writeEndElement();

	w.writeStartElement("u", "end", U1087);
	w.writeCharacters(String.valueOf(ctx.getEnd()));
	w.writeEndElement();
	
	w.writeStartElement("u", "id", U1087);
	w.writeCharacters(id);
	w.writeEndElement();


	String sn = ctx.getGenotypes().stream().filter(G->G.isHet()).map(G->G.getSampleName()).findFirst().orElse(null);
	if(sn!=null) {
		w.writeStartElement("u", "het-sample", U1087);
		w.writeCharacters(sn);
		w.writeEndElement();
		}
	sn = ctx.getGenotypes().stream().filter(G->G.isHomVar()).map(G->G.getSampleName()).findFirst().orElse(null);
	if(sn!=null) {
		w.writeStartElement("u", "homvar-sample", U1087);
		w.writeCharacters(sn);
		w.writeEndElement();
		}
	sn = ctx.getGenotypes().stream().filter(G->G.isHomRef()).map(G->G.getSampleName()).findFirst().orElse(null);
	if(sn!=null) {
		w.writeStartElement("u", "homref-sample", U1087);
		w.writeCharacters(sn);
		w.writeEndElement();
		}

	if(!ctx.isFiltered()) {
		w.writeStartElement("u", "filter", U1087);
		w.writeCharacters("PASS");
		w.writeEndElement();
		}
	else {
		for(final String f: ctx.getFilters() ) {
			w.writeStartElement("u", "filter", U1087);
			w.writeCharacters(f);
			w.writeEndElement();
			}
		}

	
	w.writeStartElement("u", "svtype", U1087);
	w.writeCharacters(svType);
	w.writeEndElement();

	w.writeStartElement("u", "svlen", U1087);
	w.writeCharacters(String.valueOf(ctx.getAttributeAsInt("SVLEN",ctx.getLengthOnReference())));
	w.writeEndElement();

	w.writeStartElement("u", "n-samples", U1087);
	w.writeCharacters(String.valueOf(ctx.getNSamples()));
	w.writeEndElement();

	w.writeStartElement("u", "het-homvar", U1087);
	w.writeCharacters(String.valueOf(ctx.getGenotypes().stream().filter(G->hasCNV(G)).count()));
	w.writeEndElement();
	
	w.writeStartElement("u", "date", U1087);
	w.writeCharacters(java.time.LocalDate.now().toString());
	w.writeEndElement();



	w.writeEndElement();
	w.writeEndElement();
	w.close();
	return sw.toString();
	}

private void instanceMain(final String args[]) {
	try {
		String vcf="";
		final Set<String> excludeIds = new HashSet<>();
		int optind=0;
		while(optind < args.length) {
			if(args[optind].equals("--vcf") && optind+1< args.length) {
				optind++;
				vcf = args[optind];
				}
			else if(args[optind].equals("--excludeids") && optind+1< args.length) {
				optind++;
				excludeIds.addAll(Files.readAllLines(Paths.get(args[optind]),java.nio.charset.Charset.defaultCharset()));
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
		out.print("\t");
		out.print("rdf");
		out.print("\t");
		out.print("svtype");
		out.println();
		try(VCFIterator r = new VCFIteratorBuilder().open(System.in)) {
			long variant_id =0L;
			final VCFHeader h = r.getHeader();
			while(r.hasNext()) {
				final VariantContext ctx = r.next();
				final String id = getVariantId(ctx);
				if(excludeIds.contains(id)) continue;
				int len = ctx.getLengthOnReference();
				if(len < minLenOnReference || len > maxLenOnReference) continue;
				final String svType = ctx.getAttributeAsString("SVTYPE","");
				if(svType==null || svType.isEmpty()) continue;
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
				sb.append("Variant location is <b>"+ctx.getContig()+":"+ctx.getStart()+"-"+ctx.getEnd()+"</b>. ");
				sb.append("SVTYPE is <b>" + svType + "</b>. ");
				sb.append("SVLEN is <b>"+ len +"</b>. ");
				sb.append("Some Samples carrying an ALT include:").append(affected.stream().limit(20).collect(Collectors.joining(" "))).append(".<br/>");

				Collections.shuffle(unaffected);
				if(unaffected.size()  > max_controls)  unaffected = unaffected.subList(0,max_controls);
				variant_id++;
				out.print(this.prefix + ctx.getContig()+"_"+ctx.getStart()+"_"+ctx.getEnd()+"_"+svType+"_"+len+"_id"+variant_id+".");
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
				out.print("\t");
				out.print(rdf(ctx));
				out.print("\t");
				out.print(svType);
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

process MERGE_KNOWN {
input:
	val(meta)
	val(L)
output:
	path("merged.bed.gz"),emit:bed
	path("merged.bed.gz.tbi")
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("htslib")}
set -o pipefail

gunzip -c ${L.join(" ")} |\
	grep -v "^#" |\
	cut -f 1-4 |\
	awk -F '\t' '(\$2 != \$3)' |\
	sort -t '\t' -T . -k1,1 -k2,2n -k3,3n --unique |\
	uniq > merged.bed

bgzip merged.bed
tabix -p bed merged.bed.gz

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">merge known files</entry>
        <entry key="number.of.files">${L.size()}</entry>
</properties>
EOF
"""
}

process SPLIT_VARIANTS {
executor "local"
afterScript "rm -rf TMP"
input:
      	val(meta)
	val(vcf)
	path(excludeids)
	val(minikit)
output:
       	path("${meta.prefix?:""}variants.tsv"),emit:output
        path("version.xml"),emit:version
script:
	def extra_filter = meta.extra_vcf_filter?:""
"""
hostname 1>&2
${moduleLoad("bcftools jvarkit")}
set -o pipefail

bcftools view "${vcf}" |\
	${isBlank(extra_filter)?"":"${extra_filter} |"} \
	java -cp  \${JVARKIT_DIST}/coverageplotter.jar:${minikit} Minikit --vcf "${vcf}" \
		${excludeids.name.equals("NO_FILE")?"":"--excludeids \"${excludeids}\""} > "${meta.prefix?:""}variants.tsv"



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
afterScript "rm -rf TMP"
memory "10g"
input:
	val(meta)
	val(reference)
	val(known)
	val(gff3)
	val(row)
output:
	path("${row.prefix}out.html"),emit:output
	path("version.xml"),emit:version
script:
	def num_cases = row.max_cases?:1000000
	def num_controls = row.max_controls?:10
	def extend = row.extend?:"5.0"
	def mapq = row.mapq?:30
"""
hostname 1>&2
${moduleLoad("jvarkit samtools")}
set -o pipefail
mkdir TMP

sort -T TMP -t '\t' -k1,1 "${row.bams}" > TMP/samples.bams.tsv

echo "${row.cases}" | tr "," "\\n" | sort | uniq > TMP/cases.txt
echo "${row.controls}" | tr "," "\\n" | sort | uniq > TMP/controls.txt


join -t '\t' -1 1 -2 1 -o "2.1" TMP/cases.txt TMP/samples.bams.tsv | shuf | head -n ${num_cases} | sort | uniq > TMP/cases2.txt
mv TMP/cases2.txt TMP/cases.txt

join -t '\t' -1 1 -2 1 -o "2.1" TMP/controls.txt TMP/samples.bams.tsv | shuf | head -n ${num_controls} | sort | uniq > TMP/controls2.txt
mv TMP/controls2.txt TMP/controls.txt


join -t '\t' -1 1 -2 1 -o "2.2" TMP/cases.txt TMP/samples.bams.tsv | sort | uniq > TMP/cases.bams.list
join -t '\t' -1 1 -2 1 -o "2.2" TMP/controls.txt TMP/samples.bams.tsv | sort | uniq > TMP/controls.bams.list


cat TMP/cases.bams.list TMP/controls.bams.list  | sort | uniq > TMP/all.bams.list

## http://www.paletton.com/#

cat << "EOF" > TMP/jeter.awk
function fr(x,y,z){ return int(y+x*(z-y));}
{SN[\$1]=((NR)*1.0);}
END{for(S in SN){f=(SN[S]/NR);printf("%s\tdisplay:block;stroke-opacity:0.8;stroke:rgb(%d,%d,%d);fill:none;\\n",S,fr(f,189,151),fr(f,0,11),fr(f,57,53));}}
EOF

awk -F '\t' -v C=red  -f TMP/jeter.awk TMP/cases.txt    >  TMP/style.css

## http://www.paletton.com/#

cat << "EOF" > TMP/jeter.awk
function fr(x,y,z){ return int(y+x*(z-y));}
{SN[\$1]=((NR)*1.0);}
END{for(S in SN){f=(SN[S]/NR);printf("%s\tdisplay:block;stroke-opacity:0.7;stroke:rgb(%d,%d,%d);fill:none;\\n",S,fr(f,46,22),fr(f,66,41),fr(f,114,85));}}
EOF

awk -F '\t' -v C=blue -f TMP/jeter.awk TMP/controls.txt >> TMP/style.css

java  -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP  -jar ${JVARKIT_DIST}/coverageplotter.jar \
	-R "${reference}" \
	--mapq "${mapq}" \
	--known '${known}' --ignore-known-containing \
	--gff3 '${gff3}' \
	--css TMP/style.css \
	--extend "${extend}" \
	--region "${row.interval}" \
	TMP/all.bams.list > TMP/jeter.html

cat << EOF > TMP/jeter.xsl
<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:svg="http://www.w3.org/2000/svg" version="1.0">
  <xsl:output method="xml" omit-xml-declaration = "yes"/>
  <xsl:template match="svg:metadata">
  <svg:metadata>
    <xsl:apply-templates select="@*|node()" /> 
    ${row.rdf}
  </svg:metadata>
  </xsl:template>

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
	<entry key="mapq">${mapq}</entry>
	<entry key="extend">${extend}</entry>
	<entry key="coverageplotter.version">\$(java -jar ${JVARKIT_DIST}/coverageplotter.jar --version)</entry>
</properties>
EOF
"""
}

process PLOT_INV {
tag "${row.prefix} ${row.interval}"
afterScript "rm -rf TMP"
memory "3g"
input:
	val(meta)
	val(reference)
	val(known)
	val(gff3)
	val(row)
output:
	path("${row.prefix}out.html"),emit:output
	path("version.xml"),emit:version
script:
	def mapq = row.mapq?:30
	def colon = row.interval.indexOf(":")
	def hyphen = row.interval.indexOf("-",colon+1)
	def contig = row.interval.substring(0,colon)
	def start = (row.interval.substring(colon+1,hyphen) as int)
	def end = (row.interval.substring(hyphen+1) as int)
	def extend = 100
	def region = end-start < extend*4 ? " --region \"${row.interval}\" ":" --region \"${contig}:${Math.max(1,start - extend)}-${start + extend}\" --region \"${contig}:${Math.max(1,end - extend)}-${end + extend}\" "

"""
hostname 1>&2
${moduleLoad("jvarkit samtools")}
set -o pipefail
mkdir TMP

sort -T TMP -t '\t' -k1,1 "${row.bams}" > TMP/samples.bams.tsv

echo "${row.cases}" | tr "," "\\n" | sort | uniq | shuf | head -n1 > TMP/cases.txt

join -t '\t' -1 1 -2 1 -o "2.2" TMP/cases.txt TMP/samples.bams.tsv | sort | uniq > TMP/cases.bams.list

test -s TMP/cases.bams.list

java  -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP  -jar ${JVARKIT_DIST}/sv2svg.jar \
	-R "${reference}" \
	--mapq "${mapq}" \
	--depth -1 \
	--mismatch \
	${region} \
	TMP/cases.bams.list > TMP/jeter.html

cat << EOF > TMP/jeter.xsl
<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:svg="http://www.w3.org/2000/svg" version="1.0">
  <xsl:output method="xml" omit-xml-declaration = "yes"/>
  <xsl:template match="svg:metadata">
  <svg:metadata>
    <xsl:apply-templates select="@*|node()" /> 
    ${row.rdf}
  </svg:metadata>
  </xsl:template>

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
	<entry key="description">plot INV</entry>
	<entry key="interval">${row.interval}</entry>
	<entry key="mapq">${mapq}</entry>
	<entry key="extend">${extend}</entry>
	<entry key="sv2svg.version">\$(java -jar ${JVARKIT_DIST}/sv2svg.jar --version)</entry>
</properties>
EOF
"""
}


process MAKEGIF {
tag "N=${L.size()}"
input:
	val(meta)
	val(L)
output:
	path("${meta.prefix?:""}coverage.gif"),emit:output
	path("version.xml"),emit:version
when:
	false
script:
"""
hostname 1>&2
${moduleLoad("ImageMagick")}

mkdir TMP
cat << EOF > TMP/jeter.list
${L.join("\n")}
EOF

i=1
cat TMP/jeter.list | while read F
do
	xmllint --xpath "//*[local-name()='svg']"   "\$F"  > TMP/jeter.svg
	convert TMP/jeter.svg TMP/jeter\${i}.png
	rm TMP/jeter.svg
	i=\$((i+1))
done

convert -resize '50%' -delay 1 -loop 0 TMP/jeter*.png "${meta.prefix?:""}coverage.gif"

rm -f TMP/*.png
###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">Convert images to GIF</entry>
</properties>
EOF
"""
}

process MERGE_PLOTS {
tag "N=${L.size()}"
input:
	val(meta)
	path(table)
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
ln -s "${table}" "${prefix}/"

${copy}

cat << EOF > ${prefix}/index.html
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
<meta http-equiv="author" content="Pierre Lindenbaum Phd ">
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

