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
nextflow.enable.dsl=2

def gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults().putReference()


gazoduc.build("bams", "NO_FILE").
	desc("File containing the paths to the indexed BAM/CRAM files.").
	existingFile().
	required().
	put()

gazoduc.make("max_cases",2).
	description("max number of cases per plot").
	setInt().
	put()

gazoduc.make("max_controls",2).
	description("max number of controls per plot").
	setInt().
	put()

gazoduc.make("vcf","NO_FILE").
	description("input vcf").
	existingFile().
	required().
	put()

gazoduc.make("extraCmdWallyRegion", "--clip").
	description("extra arguments to wally regions").
	put()

gazoduc.make("mapq", 1).
	description("min mapping quality").
	put()


include {WALLY_DOWNLOAD_01} from '../../modules/wally/wally.download.01.nf'
include {runOnComplete} from '../../modules/utils/functions.nf'
include {SAMTOOLS_SAMPLES01} from '../../modules/samtools/samtools.samples.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {moduleLoad;isBlank;isHg38;isHg19} from '../../modules/utils/functions.nf'
include {DOWNLOAD_GNOMAD_SV_01} from '../../subworkflows/gnomad/download_gnomad_sv.01.nf'
include {DOWNLOAD_DGV_01} from '../../modules/dgv/download.dgv.01.nf'
include {SIMPLE_ZIP_01} from '../../modules/utils/zip.simple.01.nf'
include {SIMPLE_PUBLISH_01} from '../../modules/utils/publish.simple.01.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'

if( params.help ) {
    gazoduc.usage().
	name("Wally").
	desc("Plot coverage for a set of bams using wally").
	print();
    exit 0
} else {
   gazoduc.validate();
}



workflow {
	ch = WALLY_REGION_01(params, params.reference, file(params.vcf), params.bams, file("NO_FILE"), file("NO_FILE") )
	html = VERSION_TO_HTML(params,ch.version)
	SIMPLE_PUBLISH_01(params, Channel.empty().mix(html.html).mix(ch.version).mix(ch.zip).collect())
	}

runOnComplete(workflow);

workflow WALLY_REGION_01 {
    take:
	    meta
	    reference
	    vcf
	    bams
	    bed
	    excludeids
    main:
		version_ch = Channel.empty()


		ch1_ch = SAMTOOLS_SAMPLES01([:],reference,bams)
		version_ch = version_ch.mix(ch1_ch.version)

	        wally_ch = WALLY_DOWNLOAD_01(meta)
                version_ch = version_ch.mix(wally_ch.version)

		compile_ch = COMPILE_VCF_PARSER(meta,reference)
		version_ch = version_ch.mix(compile_ch.version)


		merge_ch = Channel.empty()
		gnomad_ch = DOWNLOAD_GNOMAD_SV_01(meta,reference)
		version_ch = version_ch.mix(gnomad_ch.version)
		merge_ch = merge_ch.mix(gnomad_ch.bed)

		dgv_ch = DOWNLOAD_DGV_01(meta,reference)
		version_ch = version_ch.mix(dgv_ch.version)
		merge_ch = merge_ch.mix(dgv_ch.bed)

		known_ch = MERGE_KNOWN(meta,merge_ch.collect())
		version_ch = version_ch.mix(known_ch.version)


	        splitctx_ch = SPLIT_VARIANTS(meta,vcf,excludeids, compile_ch.jar)		
		version_ch = version_ch.mix(splitctx_ch.version)

		ch2_ch = splitctx_ch.output.splitCsv(header:true,sep:'\t').
			combine(ch1_ch.output).
			map{T->T[0].plus([
				"bams":T[1],
				"max_cases":(meta.max_cases?:1000000),
				"max_controls":(meta.max_controls?:10)
				])}

	        plot_ch = PLOT_WALLY(meta, reference, wally_ch.executable, known_ch.bed, ch2_ch)
		version_ch = version_ch.mix(plot_ch.version)

		zip_ch = SIMPLE_ZIP_01(params,plot_ch.output.collect())
		version_ch = version_ch.mix(zip_ch.version)


		version_ch = MERGE_VERSION(meta, "Wally", "Wally", version_ch.collect())
    emit:
	    zip = zip_ch.zip
	    version = version_ch
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

private final int minLenOnReference = ${meta.minCnvLength?:"1"};
private final int maxLenOnReference = ${meta.maxCnvLength?:"250_000_000"};
private final String prefix="${meta.prefix?:""}";
private final int max_controls =  ${meta.max_controls?:"50"};
private final int large_length =  5_000;
private final int image_size = 1024;

private boolean hasCNV(final Genotype g) {	
	return g.isHet() || g.isHomVar();
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
		out.print("title");
		out.print("\t");
		out.print("interval");
		out.print("\t");
		out.print("command");
		out.print("\t");
		out.print("vcf");
		out.print("\t");
		out.print("cases");
		out.print("\t");
		out.print("controls");
		out.print("\t");
		out.print("svtype");
		out.println();
		try(VCFIterator r = new VCFIteratorBuilder().open(System.in)) {
			long variant_id =0L;
			final VCFHeader h = r.getHeader();
			final SAMSequenceDictionary  dict = h.getSequenceDictionary();
			if(dict==null || dict.isEmpty()) {
				System.err.println("Missing dict.");
				System.exit(-1);
				}
			while(r.hasNext()) {
				final VariantContext ctx = r.next();
				final SAMSequenceRecord ssr = dict.getSequence(ctx.getContig());
				if(ssr==null) continue;
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
				

				Collections.shuffle(unaffected);
				if(unaffected.size()  > max_controls)  unaffected = unaffected.subList(0,max_controls);
				variant_id++;
				out.print(this.prefix + ctx.getContig()+"_"+ctx.getStart()+"_"+ctx.getEnd()+"_"+svType+"_" + len);
				out.print("\t");
				out.print(ctx.getContig()+":"+ctx.getStart()+"-"+ctx.getEnd());
				out.print("\t");
				if(len <= large_length*2 ) {
					final int xstart = Math.max(1, ctx.getStart() - len);
					final int xend = Math.min(ssr.getSequenceLength() , ctx.getEnd() + len);
					out.print("--region "+ctx.getContig()+":"+xstart+"-"+xend);
					out.print(" --width "+(image_size));
					}
				else
					{
					final int xstart1 = Math.max(1, ctx.getStart() - large_length);
					final int xend1 = Math.min(ssr.getSequenceLength() , ctx.getStart() + large_length);

					final int xstart2 = Math.max(1, ctx.getEnd() - large_length);
					final int xend2 = Math.min(ssr.getSequenceLength() , ctx.getEnd() + large_length);

					out.print(" --split 2 --region ");
					out.print(ctx.getContig()+":"+xstart1+"-"+xend1);
					out.print(",");
					out.print(ctx.getContig()+":"+xstart2+"-"+xend2);
					out.print(" --width "+(image_size*2));				
					}
				out.print(" --height "+ ((affected.size()+unaffected.size())*image_size)+" ");
				out.print("\t");
				out.print(vcf);
				out.print("\t");
				out.print(String.join(",",affected));
				out.print("\t");
				out.print(String.join(",",unaffected));
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


process PLOT_WALLY {
tag "${row.interval}"
afterScript "rm -rf TMP"
memory "10g"
input:
	val(meta)
	val(reference)
        path(wally)
	val(known)
	val(row)
output:
	path("${row.title}.png"),emit:output
	path("version.xml"),emit:version
script:
	def num_cases = row.max_cases?:1000000
	def num_controls = row.max_controls?:10
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


${wally.toRealPath()} region ${meta.extraCmdWallyRegion} \
	--genome "${reference}"  \
	--bed "${known}" \
	--map-qual ${meta.mapq} \
	${row.command} \
	`cat TMP/all.bams.list`

mv -v *.png "${row.title}.png"

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">plot CNV</entry>
	<entry key="interval">${row.interval}</entry>
	<entry key="mapq">${mapq}</entry>
</properties>
EOF
"""
}

