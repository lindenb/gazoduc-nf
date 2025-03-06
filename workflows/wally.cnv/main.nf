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

// use conda include {WALLY_DOWNLOAD_01} from '../../modules/wally/wally.download.01.nf'
include {runOnComplete} from '../../modules/utils/functions.nf'
include {DOWNLOAD_GNOMAD_SV_01} from '../../subworkflows/gnomad/download_gnomad_sv.01.nf'
include {DOWNLOAD_DGV_01} from '../../modules/dgv/download.dgv.01.nf'

if( params.help ) {
    exit 0
}




workflow {
	genome = Channel.of(file(params.fasta),file(params.fai),file(params.dict)).collect()
	ch = WALLY_REGION_01(genome, file(params.vcf), file(params.samplesheet) )
	}

runOnComplete(workflow);

workflow WALLY_REGION_01 {
    take:
	    genome
	    vcf
	    samplesheet
    main:

		ch1_ch = Channel.fromPath(samplesheet).splitCsv(header:true,sep:'\t')

		merge_ch = Channel.empty()

		compile_ch = COMPILE_VCF_PARSER(genome)
		gnomad_ch = DOWNLOAD_GNOMAD_SV_01(genome)
		merge_ch = merge_ch.mix(gnomad_ch.bed)

		dgv_ch = DOWNLOAD_DGV_01(genome)
		merge_ch = merge_ch.mix(dgv_ch.bed)

		known_ch = MERGE_KNOWN(merge_ch.collect())

	        splitctx_ch = SPLIT_VARIANTS(vcf,excludeids, compile_ch.jar)		

		ch2_ch = splitctx_ch.output.splitCsv(header:true,sep:'\t').
			combine(ch1_ch.rows.filter{T->T.genomeId.equals(genomeId)}).
			map{T->T[0].plus([
				"bams":T[1],
				"max_cases":(params.max_cases?:1000000),
				"max_controls":(params.max_controls?:10)
				])}

	        plot_ch = PLOT_WALLY(genome, wally_ch.executable, known_ch.bed, ch2_ch)

    emit:
	   	output = plot_ch.output
    }


process COMPILE_VCF_PARSER {
executor "local"
afterScript "rm -rf TMP"
input:
	val(genome)
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

private final int minLenOnReference = ${params.minCnvLength?:"1"};
private final int maxLenOnReference = ${params.maxCnvLength?:"250_000_000"};
private final String prefix="${params.prefix?:""}";
private final int max_controls =  ${params.max_controls?:"50"};
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
       	path("${params.prefix?:""}variants.tsv"),emit:output
script:
	def extra_filter = params.extra_vcf_filter?:""
"""
hostname 1>&2
${moduleLoad("bcftools jvarkit")}
set -o pipefail

bcftools view "${vcf}" |\
	${isBlank(extra_filter)?"":"${extra_filter} |"} \
	java -cp  \${JVARKIT_DIST}/coverageplotter.jar:${minikit} Minikit --vcf "${vcf}" \
		${excludeids.name.equals("NO_FILE")?"":"--excludeids \"${excludeids}\""} > "${params.prefix?:""}variants.tsv"

"""
}

process MERGE_KNOWN {
input:
	path("INPUT/*")
output:
	path("merged.bed.gz"),emit:bed
	path("merged.bed.gz.tbi")
script:
"""
hostname 1>&2
set -o pipefail

gunzip -c INPUT/*.gz |\
	grep -v "^#" |\
	cut -f 1-4 |\
	awk -F '\t' '(\$2 != \$3)' |\
	sort -t '\t' -T . -k1,1 -k2,2n -k3,3n --unique |\
	uniq > merged.bed

bgzip merged.bed
tabix -p bed merged.bed.gz
"""
}


process PLOT_WALLY {
tag "${row.interval}"
afterScript "rm -rf TMP"
memory "10g"
input:
	val(genome)
        path(wally)
	val(known)
	val(row)
output:
	path("${row.title}.png"),emit:output
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
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


${wally.toRealPath()} region ${params.extraCmdWallyRegion} \
	--genome "${reference}"  \
	--bed "${known}" \
	--map-qual ${params.mapq} \
	${row.command} \
	`cat TMP/all.bams.list`

mv -v *.png "${row.title}.png"
"""
}

