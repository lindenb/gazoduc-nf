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

def gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults().putReference()

gazoduc.build("vcf","NO_FILE").
	desc("vcf containing the variant").
	existingFile().
	required().
	put()

gazoduc.make("bams","NO_FILE").
        description("File containing the paths to the BAM/CRAMS files. One path per line").
        required().
        existingFile().
        put()

gazoduc.make("num_cases",-1).
        description("number of samples carrying a ALT allele to choose from the VCF. Or -1 to add all cases").
	setInt().
        put()

gazoduc.make("num_controls",-1).
        description("number of samples carrying *NO* ALT allele to choose from the VCF. Or -1 to add all controls.").
	setInt().
        put()


gazoduc.make("extraCnvTview", " --mapq 10  --cols 100 --flush").
	description("extra parameters to cnvtview").
	put();


include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {runOnComplete;moduleLoad} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {SAMTOOLS_SAMPLES_01} from '../../subworkflows/samtools/samtools.samples.01.nf'
///include {DOWNLOAD_REFGENE} from '../../modules/ucsc/download.refgene.nf'
include {SIMPLE_ZIP_01} from '../../modules/utils/zip.simple.01.nf'
include {SIMPLE_PUBLISH_01} from '../../modules/utils/publish.simple.01.nf'


if( params.help ) {
    gazoduc.usage().
	name("CNV Tview").
	desc("call jvarkit/CNVTview to a set of bam files and a vcf file").
	print();
    exit 0
} else {
   gazoduc.validate();
}


workflow {
	ch1 = CNVTVIEW(params,params.reference,file(params.bams),file(params.vcf))
	html = VERSION_TO_HTML(params,ch1.version)

	
	
	tozip = Channel.empty().mix(ch1.htmls).mix(ch1.version).mix(html.html)
	zip_ch = SIMPLE_ZIP_01([:] ,tozip.collect())

	SIMPLE_PUBLISH_01(params, Channel.empty().mix(zip_ch.zip).collect())
	}

runOnComplete(workflow)

workflow CNVTVIEW {
	take:
		meta
		reference
		bams
		vcf
	main:
		version_ch = Channel.empty()

                snbam_ch = SAMTOOLS_SAMPLES_01(meta.plus("with_header":false,"allow_multiple_references":false,"allow_duplicate_samples":false), reference, file("NO_FILE"), bams)
                version_ch = version_ch.mix(snbam_ch.version)

		compile_ch = COMPILE(meta)
		version_ch = version_ch.mix(compile_ch.version)

		//refgene_ch = DOWNLOAD_REFGENE(meta,reference)
                //version_ch = version_ch.mix(refgene_ch.version)

		prepare_ch = PREPARE_CNVTVIEW(meta, reference, compile_ch.output, snbam_ch.output, vcf )
                version_ch = version_ch.mix(prepare_ch.version)

		report_ch = INVOKE_CNVTVIEW(meta,reference, prepare_ch.output.splitCsv(header:true,sep:'\t') )
                version_ch = version_ch.mix(report_ch.version)
		
		version_ch = MERGE_VERSION(meta, "IGV report", "IGV report", version_ch.collect())

	emit:
		version = version_ch
		htmls  = report_ch.output.collect()
	}


process COMPILE {
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

cat << EOF > Minikit.java
import java.io.BufferedReader;
import java.io.PrintWriter;
import java.nio.file.*;
import java.util.*;
import java.util.function.ToIntFunction;
import java.util.logging.Logger;
import java.util.regex.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

public class Minikit {
	private static final Logger LOG = Logger.getLogger("Minikit");
	private final int MAX_CONTROLS= ${meta.num_controls};
	private final int MAX_CASES= ${meta.num_cases};
	private final Map<String,Path> sample2bam = new HashMap<>();
	
	private List<Locatable> parseBnd(final VariantContext ctx,final SAMSequenceDictionary dict) {
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
                final String contig =S.substring(x1+1,colon);
                if(dict.getSequence(contig)==null) return null;
                int pos = Integer.parseInt(S.substring(colon+1,x2));
                return new Interval(contig,pos,pos);
                }).
                filter(R->R!=null).
                collect(Collectors.toList());
			}

	
    public int instanceMainWithExit(String[] args) {
    	final Pattern tab = Pattern.compile("\\t");
        int optind=0;
        Path bams = null;
        Path outdir = null;
        SAMSequenceDictionary dict = null;
        while(optind < args.length) {
                if(args[optind].equals("--bams") && optind+1< args.length) {
                    optind++;
                    bams = Paths.get(args[optind]);
                    }
                else if(args[optind].equals("--reference") && optind+1< args.length) {
                    optind++;
                    dict = SAMSequenceDictionaryExtractor.extractDictionary(Paths.get(args[optind]));
                    }
                else if(args[optind].equals("--out") && optind+1< args.length) {
                    optind++;
                    outdir = Paths.get(args[optind]);
                    IOUtil.assertDirectoryIsWritable(outdir);
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
        if(bams==null) {
			System.err.println("undefined bams.");
			System.exit(-1);
        	}
        if(outdir==null) {
			System.err.println("undefined outdir.");
			System.exit(-1);
        	}
        if(dict==null) {
			System.err.println("undefined dict.");
			System.exit(-1);
        	}
    	try                              
            {
    	    final int SV_FLANKING= ${meta.sv_flanking?:250};
            try(BufferedReader in=Files.newBufferedReader(bams)) {
                    in.lines().
                    	filter(T->!T.startsWith("#")).
                    	map(T->tab.split(T)).
                    	forEach(T->sample2bam.put(T[0],Paths.get(T[2])));
                    }
            final ToIntFunction<Genotype> countALT = G->(int)(G.getAlleles().stream().filter(A->!(A.isReference() || A.isNoCall())).count());
            final Comparator<Genotype> gtSorter = (A,B)->{
            	if(A.hasGQ() && B.hasGQ()) {
            		return Integer.compare(B.getGQ(), A.getGQ());
            		}
            	if(A.hasDP() && B.hasDP())  {
            		return Integer.compare(B.getDP(), A.getDP());
            		}
            	return A.getSampleName().compareTo(B.getSampleName());
            	};
            	
            final List<String> output_keys = Arrays.asList("title","interval","bams","top","header");
            System.out.println(String.join("\\t", output_keys));
            try(VCFIterator r = new VCFIteratorBuilder().open(System.in)) {
                final VCFHeader header = r.getHeader();
                while(r.hasNext()) {
                        final VariantContext ctx = r.next();
                        final String svType = ctx.getAttributeAsString("SVTYPE","");
                        if(svType.equals("BND")) continue;
                        if(svType.equals("INV")) continue;
			final int chromLen = dict.getSequence(ctx.getContig()).getSequenceLength();

			final int dx5;
			final int dx3;
			if(ctx.hasAttribute("CIPOS")) {
				dx5 = ctx.getAttributeAsIntList("CIPOS",0).stream().mapToInt(V->Math.abs(V)).max().orElse(0);
			} else {
				dx5 = 0;
			}
			if(ctx.hasAttribute("CIEND")) {
				dx3 = ctx.getAttributeAsIntList("CIEND",0).stream().mapToInt(V->Math.abs(V)).max().orElse(0);
			} else {
				dx3 = 0;
			}


                        final List<Genotype> cases = ctx.getGenotypes().stream().
                            	filter(G->sample2bam.containsKey(G.getSampleName())).
                            	filter(G->countALT.applyAsInt(G)>0L).
                            	sorted(gtSorter).
                            	limit(MAX_CASES==-1?Integer.MAX_VALUE:MAX_CASES).
                            	collect(Collectors.toList());
                        if(cases.isEmpty()) continue;
                        
                        final List<Genotype> controls = ctx.getGenotypes().stream().
                        	filter(G->sample2bam.containsKey(G.getSampleName())).
                        	filter(G->countALT.applyAsInt(G)==0L).
                        	sorted(gtSorter).
                        	limit(MAX_CONTROLS==-1?Integer.MAX_VALUE:MAX_CONTROLS).
                        	collect(Collectors.toList());

                        
                        final Set<String> sample_names = Stream.concat(cases.stream(),controls.stream()).
                        		map(G->G.getSampleName()).
                        		collect(Collectors.toSet());
                        
                        final Map<String,String> properties = new HashMap<>();
			final List<String> headers = new ArrayList<>(10);

			headers.add("Interval   : " + ctx.getContig()+":"+ctx.getStart()+"-"+ctx.getEnd());
			headers.add("SVLEN      : " + ctx.getLengthOnReference());
			headers.add("SVTYPE     : " + svType);
			headers.add("CASES      : " + cases.stream().map(G->G.getSampleName()).sorted().collect(Collectors.joining(" ")));
			headers.add("N-CASES    : " + cases.size());
			headers.add("CONTROLS   : " + controls.stream().map(G->G.getSampleName()).sorted().collect(Collectors.joining(" ")));
			headers.add("N-CONTROLS : " + controls.size());

                        StringBuilder vcNameBuilder  = new StringBuilder(ctx.getContig()).append("_").append(ctx.getStart());
                        if(ctx.getEnd()!=ctx.getStart()) vcNameBuilder.append("_").append(ctx.getEnd()).append("_").append(ctx.getLengthOnReference());
                        if(!svType.isEmpty())  vcNameBuilder.append("_").append(svType);
                        final String vcName= vcNameBuilder.toString();
                        properties.put("title",vcName);


			properties.put("header",String.join("|",headers));

                        properties.put("bams",sample_names.stream().
                        	map(S->sample2bam.get(S).toString()).
                        	collect(Collectors.joining(" "))
				);
                        
			properties.put("top",
				cases.
				stream().
				map(G->G.getSampleName()).
				sorted().
				collect(Collectors.joining(","))
				);

			properties.put("interval",
                        	ctx.getContig() + ":" +
                        	Math.max(0, ctx.getStart()-(dx5+1)) + "-" +
                        	Math.min(chromLen, ctx.getEnd()+dx3)
				);
                        System.out.println(output_keys.stream().map(K->properties.getOrDefault(K,"")).collect(Collectors.joining("\t")));
                    	}
            	}
            return 0;                                                                          
            }             
    catch(Throwable err)
            {                                                                                  
            err.printStackTrace();
            return -1;
            }                            
    	}                                  


public static void main(String[] args)
    {
    System.exit(new Minikit().instanceMainWithExit(args));
    }

}
EOF

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

process PREPARE_CNVTVIEW {
input:
	val(meta)
	val(reference)
	path(jar)
	path(snbam)
	path(vcf)
output:
	path("output.tsv"),emit:output
	path("version.xml"),emit:version
script:
"""        
hostname 1>&2     
${moduleLoad("jvarkit bcftools")}
set -o pipefail
mkdir -p OUT

bcftools view "${vcf}" |\
	java -cp ${jar}:\${JVARKIT_DIST}/vcffilterjdk.jar Minikit --reference "${reference}" --bams "${snbam}" --out "\${PWD}/OUT" > output.tsv

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">prepare data for CNVTVIEW</entry>
        <entry key="versions"></entry>
</properties>
EOF
"""
}

process INVOKE_CNVTVIEW {
tag "${row.title}"
memory "5g"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
	val(row)
output:
	path("${row.title}.txt"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("jvarkit")}
mkdir -p TMP

echo '${row.bams}' | tr " " "\\n" | sort | uniq > TMP/jeter.list


echo '${row.header}' | tr '|' '\\n' > TMP/jeter.txt
echo >> TMP/jeter.txt

java -jar -Xmx${task.memory.giga}G  -Djava.io.tmpdir=TMP \${JVARKIT_DIST}/cnvtview.jar \
	-R "${reference}" \
	--interval "${row.interval}" \
	${meta.extraCnvTview} \
	${row.top.isEmpty() || meta.extraCnvTview.contains("--flush")?"":row.top.split(",").collect{x->"--top "+x}.join(" ")} \
	TMP/jeter.list >> TMP/jeter.txt


mv -v "TMP/jeter.txt" "${row.title}.txt" 

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">generate CNVTview reports</entry>
        <entry key="versions"></entry>
</properties>
EOF
"""
}
