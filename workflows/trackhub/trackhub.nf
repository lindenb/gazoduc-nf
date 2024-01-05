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
nextflow.enable.dsl=2

def gazo = gazoduc.Gazoduc.getInstance(params).putDefaults().putGenomes()

include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {escapeXml;runOnComplete;moduleLoad} from '../../modules/utils/functions.nf'
include {SIMPLE_PUBLISH_01} from '../../modules/utils/publish.simple.01.nf'


gazo.build("vcfs","NO_FILE").
	desc(gazoduc.Gazoduc.DESC_VCF_LIST).
	existingFile().
	required().
	put()

gazo.build("with_bnd",true).
	desc("include variants with INFO/SVTYPE=BND").
	setBoolean().
	put()

gazo.build("min_sv_len",0).
	desc("min SV length").
	setInteger().
	put()


gazo.build("max_sv_len",10_000_000).
	desc("max SV length").
	setInteger().
	put()

gazo.build("description","").
	desc("track description").
	required().
	notEmpty().
	put()



if( params.help ) {
    gazo.usage().
	name("trackhub").
	desc("Build a UCSC track hub form a set of SV/CNV VCFs files.").
	print();
    exit 0
} else {
   gazo.validate();
}

workflow {
	ch1 = TRACK_HUB(params, file(params.vcfs))
	html = VERSION_TO_HTML(params,ch1.version)

	pub_ch = Channel.empty().mix(ch1.tar).mix(ch1.version).mix(html.html)
	SIMPLE_PUBLISH_01(params, pub_ch.collect())
	}

runOnComplete(workflow)

workflow TRACK_HUB {
	take:
		meta
		vcfs
	main:
		version_ch = Channel.empty()

		c1_ch = COMPILE_JAVA(meta)
                version_ch = version_ch.mix(c1_ch.version)
		
		ci_ch = DOWNLOAD_CHROM_INFO(meta)
		version_ch = version_ch.mix(ci_ch.version)

		scan_ch = SCAN_VCF( meta, c1_ch.output, Channel.fromPath(vcfs).splitText().map{file(it.trim())} )
		version_ch = version_ch.mix(scan_ch.version)

		types_ch = ALL_TYPES(meta, scan_ch.types.collect())
		version_ch = version_ch.mix(types_ch.version)

		each_type = types_ch.output.splitText().map{it.trim()}


		all_outputs = Channel.empty()
		bed_ch = MERGE_BED(meta, ci_ch.output , vcfs,  each_type.combine(scan_ch.bed).groupTuple() )
		version_ch = version_ch.mix(bed_ch.version)
		all_outputs = all_outputs.mix(bed_ch.output)
	
		if ( meta.with_bnd as boolean ) {
			inter_ch = MERGE_BED_PE(meta, ci_ch.output , vcfs, scan_ch.bedpe.collect())
			version_ch = version_ch.mix(inter_ch.version)
			all_outputs = all_outputs.mix(inter_ch.output)
			}

		track_ch = MAKE_TRACK(meta, all_outputs.collect() )
		version_ch = version_ch.mix(track_ch.version)

		version_ch = MERGE_VERSION(meta, "trackhub", "Trackhub", version_ch.collect())
	emit:
		version = version_ch
		tar = track_ch.tar
	}



process COMPILE_JAVA {
	afterScript "rm -rf TMP"
	input:
		val(meta)	
	output:
		path("minikit.jar"),emit:output
		path("version.xml"),emit:version
	script:
	"""
	hostname 1>&2
	${moduleLoad("jvarkit")}

cat << __EOF__ > Minikit.java
import java.nio.file.*;
import java.io.*;
import java.util.*;
import java.util.stream.*;

import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.*;
import htsjdk.variant.vcf.*;


public class Minikit {
private static final boolean WITH_BND = ${meta.with_bnd as boolean};
private final Map<String,String> convertHash = new HashMap<>();
private static final int MIN_SV_LEN = ${meta.min_sv_len};
private static final int MAX_SV_LEN = ${meta.max_sv_len};

 
private String convert(final String s) {
	return convertHash.get(s);
	}
private String type2color(String s) {
	if(s.equals("DEL")) return "0,0,255";
	if(s.equals("INS")) return "0,255,0";
	if(s.equals("INV")) return "255,0,0";
	if(s.equals("BND")) return "255,255,0";
	if(s.contains("DUP")) return "255,0,255";
	return "0,0,0";
	}
// A]hs37d5:12060965]]
private List<Locatable> parseBnd(final VariantContext ctx) {
	if(!WITH_BND) return Collections.emptyList();
	return ctx.getAlternateAlleles().stream().
		filter(A->A.isSymbolic()).
		map(A->A.getDisplayString()).
		flatMap(S->Arrays.stream(S.split("[,]"))).
		map(S->{
		int x1 = S.indexOf("[");
		if(x1==-1) x1= S.indexOf("]");
		if(x1==-1) return null;
		int x2 =  S.indexOf("[",x1+1);
		if(x2==-1) x2 =S.indexOf("]",x1+1);
		if(x2==-1) return null;
		int colon = S.indexOf(":",x1+1);
		if(colon==-1 || colon>=x2) return null;
		final String contig = convert(S.substring(x1+1,colon));
		if(contig==null) return null;
		int pos = Integer.parseInt(S.substring(colon+1,x2));
		return new Interval(contig,pos,pos);
		}).
		filter(R->R!=null).
		collect(Collectors.toList());
	}

private int doWork(final List<String> args) {
	if(args.size()!=3) {
		System.err.println("usage vcfname out.bed out.bedpe");
		System.exit(-1);
		}
	try 
		{
		
		String sampleName = null;
		File out1 = null;
		File out2 = null;
		int optind=0;
		for(int i=1;i<=22;i++) {
			convertHash.put(""+i,"chr"+i);
			convertHash.put("chr"+i,"chr"+i);
			}

		convertHash.put("X","chrX");
		convertHash.put("chrX","chrX");
		convertHash.put("Y","chrY");
		convertHash.put("chrY","chrY");
		convertHash.put("MT","chrM");
		convertHash.put("M","chrM");

		String vcfName= new File(args.get(0)).getName();
		if(vcfName.endsWith(".vcf.gz")) vcfName=vcfName.substring(0, vcfName.length() -7);
		if(vcfName.endsWith(".vcf")) vcfName=vcfName.substring(0, vcfName.length() -4);
		if(vcfName.endsWith(".bcf")) vcfName=vcfName.substring(0, vcfName.length() -4);

		try(PrintWriter pw1 = new PrintWriter(args.get(1)); PrintWriter pw2 = new PrintWriter(args.get(2))){
			try(final VCFIterator r = new VCFIteratorBuilder().open(System.in)) {
				final List<String> samples = r.getHeader().getSampleNamesInOrder();
				final SAMSequenceDictionary dict = r.getHeader().getSequenceDictionary();
				
				while(r.hasNext()) {
					final VariantContext ctx = r.next();
					final SAMSequenceRecord ssr = (dict==null?null:dict.getSequence(ctx.getContig()));
					if(ssr==null) continue;
					final int cipos = ctx.getAttributeAsIntList("CIPOS",0).stream().mapToInt(V->V.intValue()).min().orElse(0);
					final int ciend = ctx.getAttributeAsIntList("CIEND",0).stream().mapToInt(V->V.intValue()).max().orElse(0);

					int start0 = Math.max(0,(ctx.getStart()-1) + cipos);
					int end0 = ctx.getEnd() + ciend;

					if(start0>=end0 || end0>ssr.getSequenceLength()) {
						start0 = ctx.getStart()-1;
						end0 = ctx.getEnd();
						}

					final String type= ctx.getAttributeAsString("SVTYPE",null);
					if(type==null || type.isEmpty()) continue;
					String contig =  convert(ctx.getContig());
					if(contig==null || contig.isEmpty()) continue;

					for(Genotype gt : ctx.getGenotypes()) {
						if(gt.isNoCall() || gt.isHomRef()) continue;
						if(type.equals("BND")) {
						for(Locatable bnd : parseBnd(ctx)) {
							final String contig2 =  convert(bnd.getContig());
							if(contig2==null || contig2.isEmpty()) continue;
							pw2.print(contig);
							pw2.print("\t");
							pw2.print(start0);
							pw2.print("\t");
							pw2.print(end0);
							pw2.print("\t");
							pw2.print(gt.getSampleName());
							pw2.print("\t");
							pw2.print(ctx.isFiltered()?"100":"900");
							pw2.print("\t");
							pw2.print(ctx.isFiltered()?"1":"2");//value
							pw2.print("\t");
							pw2.print(gt.getSampleName());//exp
							pw2.print("\t");
							pw2.print(contig.equals("contig2")?"0,0,255":"255,0,0");//color
							pw2.print("\t");
							pw2.print(contig);//source chrom
							pw2.print("\t");
							pw2.print(ctx.getStart()-1);//source start
							pw2.print("\t");
							pw2.print(ctx.getEnd());//source end
							pw2.print("\t");
							pw2.print(".");//source name
							pw2.print("\t");
							pw2.print(".");//source strand
							pw2.print("\t");
							pw2.print(contig2);//target chrom
							pw2.print("\t");
							pw2.print(bnd.getStart()-1); //target start
							pw2.print("\t");
							pw2.print(bnd.getEnd());//target end
							pw2.print("\t");
							pw2.print(".");//target name
							pw2.print("\t");
							pw2.print(".");//target strand
							pw2.print("\t");
							pw2.print(ctx.isFiltered()?String.join("_",ctx.getFilters()):".");//filters
							pw2.print("\t");
							pw2.print(gt.getType().name());//gt type
							pw2.print("\t");
							pw2.print(vcfName);//vcf Name
							pw2.println();

							}
						} else
						{
						final int svlen ;
						if(ctx.hasAttribute("SVLEN")) {
							svlen = Math.abs( ctx.getAttributeAsInt("SVLEN",0));
							}
						else if(ctx.hasAttribute("SVINSLEN")) {
							svlen = Math.abs( ctx.getAttributeAsInt("SVINSLEN",0));
							}
						else if(ctx.hasAttribute("SVINSSEQ")) {
							svlen = Math.abs( ctx.getAttributeAsString("SVINSSEQ","").length());
							}
						else
							{
							svlen = ctx.getLengthOnReference();
							}
						
						if( svlen < MIN_SV_LEN) continue;
						if( svlen > MAX_SV_LEN) continue;

						pw1.print(contig);
						pw1.print("\t");
						pw1.print(start0);
						pw1.print("\t");
						pw1.print(end0);
						pw1.print("\t");
						pw1.print(gt.getSampleName());
						pw1.print("\t");
						pw1.print(ctx.isFiltered()?"100":"900");
						pw1.print("\t");
						pw1.print("+");
						pw1.print("\t");
						pw1.print(ctx.getStart()-1);
						pw1.print("\t");
						pw1.print(ctx.getEnd());
						pw1.print("\t");
						pw1.print(type2color(type));
						pw1.print("\t");
						pw1.print(type);//type
						pw1.print("\t");
						pw1.print(svlen);//svlen
						pw1.print("\t");
						pw1.print(ctx.isFiltered()?String.join("_",ctx.getFilters()):".");//filters
						pw1.print("\t");
						pw1.print(gt.getType().name());//gt type
						pw1.print("\t");
						pw1.print(vcfName);//vcf Name
						pw1.println();
						}
					} // end loop genotypes
				}
			
			}
		pw1.flush();
		pw2.flush();
		}
	catch(final Throwable err) {
		err.printStackTrace();
		return -1;
		}

	
	return 0;
	}
catch(final Throwable err2) {
	err2.printStackTrace();
	return -1;	
	}
}

public static void main(final String args[]) {
	try {
		int ret= new Minikit().doWork(Arrays.asList(args));
		System.exit(ret);
		}
	catch(Throwable err) {
		err.printStackTrace();
		System.exit(-1);
		}
	}

}
__EOF__

mkdir -p TMP
javac -d TMP -cp \${JVARKIT_DIST}/vcffilterjdk.jar:. Minikit.java
jar cvf minikit.jar -C TMP .


###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">compile minikit</entry>
	<entry key="javac.version">\$(javac -version 2>&1)</entry>
</properties>
EOF
	"""
	}


process SCAN_VCF {
	tag "${vcf.name}"
	afterScript "rm -rf TMP"
	input:
		val(meta)
		path(minikit)
		path(vcf)
	output:
		path("all.bed"),emit:bed
		path("all.bedpe"),emit:bedpe
		path("types.txt"),emit:types
		path("version.xml"),emit:version
	when:
		true
	script:
	"""
	hostname 1>&2
	${moduleLoad("bcftools jvarkit")}
	mkdir -p TMP
	
	bcftools view "${vcf}" |\
		java -Djava.io.tmpdir=TMP -cp \${JVARKIT_DIST}/vcffilterjdk.jar:${minikit} Minikit "${vcf}" TMP/tmp.bed TMP/tmp.bedpe

	cut -f 10 TMP/tmp.bed | uniq | LC_ALL=C sort -T TMP | LC_ALL=C uniq > types.txt

	LC_ALL=C sort -T TMP -k1,1 -k2,2n TMP/tmp.bed | uniq > all.bed
	LC_ALL=C sort -T TMP -k1,1 -k2,2n TMP/tmp.bedpe | uniq > all.bedpe
	###############################################################################
	cat << EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">extract variants</entry>
	</properties>
	EOF
	"""
	}




process DOWNLOAD_CHROM_INFO {
executor "local"
input:
	val(meta)
output:
	path("chromInfo.txt"),emit:output
	path("version.xml"),emit:version
script:
	def build = gazo.getGenome().getDictionary().getUcscName().orElse("undefined");
	def url = "https://hgdownload.cse.ucsc.edu/goldenpath/${build}/database/chromInfo.txt.gz"
"""
wget -O "chromInfo.txt.gz" "${url}"
gunzip chromInfo.txt.gz

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">download chrominfo</entry>
	<entry key="url">${url}</entry>
</properties>
EOF
"""
}

process ALL_TYPES {
	tag "N=${L.size()}"
        afterScript "rm -rf TMP"
        input:
                val(meta)
                val(L)
	output:
		path("types.txt"),emit:output
		path("version.xml"),emit:version
	script:
	"""
	hostname 1>&2
	mkdir -p TMP

cat << EOF | tr "\\n"  "\\0" > TMP/files.0
${L.join("\n")}
EOF

	LC_ALL=C sort -T TMP  --files0-from=TMP/files.0 --merge | uniq > types.txt

	###############################################################################
	cat << EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">distinct SV types.</entry>
	</properties>
	EOF


	"""
	}

process MERGE_BED {
	tag "${type} N=${beds.size()}"
	afterScript "rm -rf TMP"
	input:
		val(meta)
		path(chrominfo)
		path(vcfs)
		tuple val(type),val(beds)
	output:
		path("output.tsv"),emit:output
		path("version.xml"),emit:version
	script:
		def filename = file("${meta.prefix?:""}sv.${type}.bb")
		def now =  new java.text.SimpleDateFormat("yyyy-MM-dd").format(new java.util.Date());

	"""
	hostname 1>&2
	${moduleLoad("ucsc")}
	set -o pipefail
	mkdir -p TMP

cat << EOF > TMP/schema.as
table sv
"Structural Variant"
    (
    string chrom;      "Chromosome"
    uint chromStart;   "Start position"
    uint chromEnd;     "End position"
    string name;       "Sample Name"
    uint score;        "Score from 0-1000."
    char[1] strand;    "+ or -"
    uint thickStart;   "VCF/POS"
    uint thickEnd;   "VCF/END"
    uint reserved;        "rgb color of item"
    string type; "Sv type=${type}"
    int svLen; "Sv Len"
    string filters; "FILTERs"
    string gtType; "Genotype Type"
    string vcfName; "VCF source"
    )
EOF

cat << EOF | tr "\\n"  "\\0" > TMP/files.0
${beds.join("\n")}
EOF

	LC_ALL=C sort -T TMP -k1,1 -k2,2n  --files0-from=TMP/files.0 --merge | awk -F '\t' '(\$10=="${type}")' |  uniq > TMP/all.bed

	bedToBigBed -as=TMP/schema.as -type=bed9+5 TMP/all.bed "${chrominfo}" "${filename.name}"



#################################################

cat << EOF > trackDb.txt
track ${filename.getBaseName()}
bigDataUrl ${filename.name}
shortLabel ${filename.getBaseName()}
longLabel ${filename.getBaseName()} ${type} (${meta.min_sv_len} < SVLEN < ${meta.max_sv_len})
type bigBed 9
itemRgb "On"

EOF

#################################################

cat << EOF > TMP/jeter.html
<html>
<body>
<h1>${filename.name}</h1>
<p>${type} (${meta.min_sv_len} < SVLEN < ${meta.max_sv_len}</p>
<p>${escapeXml(meta.description)}</p>
<h3>Schema</h3><pre>
EOF

cat TMP/schema.as >> TMP/jeter.html

cat << EOF >> TMP/jeter.html
</pre><h3>Sources</h3>
<pre>
EOF

cat '${vcfs}' >> TMP/jeter.html

cat << EOF >> TMP/jeter.html
</pre>
<hr/>
<div>Author: Pierre Lindenbaum. Generated on ${now}</div>
</body>
</html>
EOF

mv -v TMP/jeter.html "${filename.getBaseName()}.html"



	echo "\${PWD}/trackDb.txt\t\${PWD}/${filename.getBaseName()}.html\t\${PWD}/${filename.name}" > output.tsv


	###############################################################################
	cat << EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="type">${type}</entry>
		<entry key="description">create SV.as</entry>
	</properties>
	EOF
	"""
	}




process MERGE_BED_PE {
	tag "N=${beds.size()}"
	afterScript "rm -rf TMP"
	input:
		val(meta)
		path(chrominfo)
		path(vcfs)
		val(beds)
	output:
		path("output.tsv"),emit:output
		path("version.xml"),emit:version
	script:
		def filename = file("${meta.prefix?:""}interact.bb")
		def now =  new java.text.SimpleDateFormat("yyyy-MM-dd").format(new java.util.Date());
	"""
	hostname 1>&2
	${moduleLoad("ucsc/0.0.0")}
	mkdir -p TMP


cat << EOF > TMP/schema.as
table interact
"BND"
    (
    string chrom;      "Chromosome (or contig, scaffold, etc.). For interchromosomal, use 2 records"
    uint chromStart;   "Start position of lower region. For interchromosomal, set to chromStart of this region"
    uint chromEnd;     "End position of upper region. For interchromosomal, set to chromEnd of this region"
    string name;       "Name of item, for display.  Usually 'sourceName/targetName' or empty"
    uint score;        "Score from 0-1000."
    double value;      "Strength of interaction or other data value. Typically basis for score"
    string exp;        "Experiment name (metadata for filtering). Use . if not applicable"
    string color;      "Item color.  Specified as r,g,b or hexadecimal #RRGGBB or html color name, as in //www.w3.org/TR/css3-color/#html4."
    string sourceChrom;  "Chromosome of source region (directional) or lower region. For non-directional interchromosomal, chrom of this region."
    uint sourceStart;  "Start position source/lower/this region"
    uint sourceEnd;    "End position in chromosome of source/lower/this region"
    string sourceName;  "Identifier of source/lower/this region"
    string sourceStrand; "Orientation of source/lower/this region: + or -.  Use . if not applicable"
    string targetChrom; "Chromosome of target region (directional) or upper region. For non-directional interchromosomal, chrom of other region"
    uint targetStart;  "Start position in chromosome of target/upper/this region"
    uint targetEnd;    "End position in chromosome of target/upper/this region"
    string targetName; "Identifier of target/upper/this region"
    string targetStrand; "Orientation of target/upper/this region: + or -.  Use . if not applicable"
    string filters; "FILTER column"
    string gtType; "Genotype Type"
    string vcfName; "VCF source"
    )
EOF

cat << EOF | tr "\\n"  "\\0" > TMP/files.0
${beds.join("\n")}
EOF

	LC_ALL=C sort -T TMP -k1,1 -k2,2n --merge --files0-from=TMP/files.0 | uniq > TMP/all.bedpe

	bedToBigBed -as=TMP/schema.as -type=bed5+16 TMP/all.bedpe "${chrominfo}"  "${meta.prefix?:""}interact.bb"



cat << EOF > trackDb.txt

track ${filename.getBaseName()}
bigDataUrl ${filename.name}
shortLabel ${filename.getBaseName()}
longLabel ${filename.getBaseName()} BND
type bigInteract
itemRgb "On"

EOF


#################################################

cat << EOF > TMP/jeter.html
<html>
<h1>${filename.name}</h1>
<p>BND</p>
<p>${escapeXml(meta.description)}</p>
<h3>Schema</h3><pre>
EOF

cat "TMP/schema.as" >> TMP/jeter.html

cat << EOF >> TMP/jeter.html
</pre><h3>Sources</h3>
<pre>
EOF

cat '${vcfs}' >> TMP/jeter.html

cat << EOF >> TMP/jeter.html
</pre>
<hr/>
<div>Author: Pierre Lindenbaum. Generated on ${now}</div>
</html>
EOF


mv -v "TMP/jeter.html" "${filename.getBaseName()}.html"


	echo "\${PWD}/trackDb.txt\t\${PWD}/${filename.getBaseName()}.html\t\${PWD}/${filename.name}" > output.tsv


	###############################################################################
	cat << EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
	</properties>
	EOF
	"""
	}

process MAKE_TRACK {
	tag "N=${L.size()}"
	//afterScript "rm -rf TMP"
	input:
		val(meta)
		val(L)
	output:
		path("${meta.prefix?:""}hub.tar.gz"),emit:tar
		path("version.xml"),emit:version
	script:
		def build = gazo.getGenome().getDictionary().getUcscName().orElse("undefined")
		def prefix= meta.prefix?:""
		def now =  new java.text.SimpleDateFormat("yyyy-MM-dd").format(new java.util.Date());
	"""
	hostname 1>&2
	${moduleLoad("ucsc/0.0.0")}
	mkdir -p "TMP/${build}"
	set -x

	cat ${L.join(" ")} > TMP/jeter.list

#################################################

cat << EOF > TMP/genomes.txt
genome ${build}
trackDb ${build}/trackDb.txt
EOF

#################################################

cat << EOF > TMP/hub.html
<html>
<head><title>${prefix}hub</title></head>
<body>
<h1>${prefix}hub</h1>
<p>${escapeXml(meta.description)}</p>
<hr/>
<div>Author: Pierre Lindenbaum. Generated on ${now}</div>
</body>
</html>
EOF

#################################################

cat << EOF > TMP/hub.txt
hub ${prefix}hub
shortLabel ${prefix}hub
longLabel ${prefix}hub
genomesFile genomes.txt
email plindenbaum@yahoo.fr
descriptionUrl hub.html 
EOF

touch TMP/${build}/trackDb.txt

cut -f 1 TMP/jeter.list | xargs -L 50 cat  >> TMP/${build}/trackDb.txt

cut -f 2 TMP/jeter.list	| xargs -L1 -I '{}' cp -v '{}' TMP/${build}/
cut -f 3 TMP/jeter.list	| xargs -L1 -I '{}' cp -v '{}' TMP/${build}/

#################################################

find TMP -type f 1>&2

hubCheck "\${PWD}/TMP/hub.txt"

rm TMP/jeter.list
mv -v TMP "${prefix}hub"

tar cvfz "${prefix}hub.tar.gz" "${prefix}hub"

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
</properties>
EOF
"""
}
