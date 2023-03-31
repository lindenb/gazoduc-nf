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

gazoduc.make("max_cases",100_000).
	description("max number of cases per plot").
	setInt().
	put()

gazoduc.make("max_controls",10).
	description("max number of controls per plot").
	setInt().
	put()





include {WALLY_DOWNLOAD_01} from '../../modules/wally/wally.download.01.nf'
include {runOnComplete} from '../../modules/utils/functions.nf'
include {SAMTOOLS_SAMPLES01} from '../../modules/samtools/samtools.samples.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {moduleLoad;isBlank;isHg38;isHg19} from '../../modules/utils/functions.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {DOWNLOAD_GNOMAD_SV_01} from '../../modules/gnomad/download.gnomad.sv.01.nf'
include {DOWNLOAD_DGV_01} from '../../modules/dgv/download.dgv.01.nf'
include {DOWNLOAD_GFF3_01} from '../../modules/gff3/download.gff3.01.nf'


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

	//html = VERSION_TO_HTML(params,ch1.version)
	//PUBLISH(ch1.version,html.html,ch1.zip)
	}

runOnComplete(workflow);

workflow WALLY_REGION_01 {
    take:
	    meta
	    reference
	    bams
	    bed
    main:
		version_ch = Channel.empty()

		ch1_ch = SAMTOOLS_SAMPLES01([:],reference,bams)
		version_ch = version_ch.mix(ch1_ch.version)

	    wally_ch = WALLY_DOWNLOAD_01(meta)
        version_ch = version_ch.mix(ch1_ch.version)

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

        plot_ch = PLOT_WALLY(meta, reference, ch2_ch)
		version_ch = version_ch.mix(plot_ch.version)

		version_ch = MERGE_VERSION(meta, "Wally", "Wally", version_ch.collect())
    emit:
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




process PLOT_WALLY {
tag "${row.prefix} ${row.interval}"
afterScript "rm -rf TMP"
memory "10g"
input:
	val(meta)
	val(reference)
    path(wally)
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

