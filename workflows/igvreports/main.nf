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


include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {dumpParams;runOnComplete;moduleLoad} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {SAMTOOLS_SAMPLES} from '../../subworkflows/samtools/samtools.samples.03.nf'
include {DOWNLOAD_CYTOBAND} from '../../modules/ucsc/download.cytoband.nf'
include {DOWNLOAD_REFGENE} from '../../modules/ucsc/download.refgene.nf'

if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}



workflow {
	ch1 = IGVREPORTS([:],params.genomeId,file(params.bams),file(params.vcf))
	html = VERSION_TO_HTML(ch1.version)
	}

runOnComplete(workflow)

workflow IGVREPORTS {
	take:
		meta
		genomeId
		bams
		vcf
	main:
		version_ch = Channel.empty()

                snbam_ch = SAMTOOLS_SAMPLES(["with_header":false,"allow_multiple_references":false,"allow_duplicate_samples":false], bams)
                version_ch = version_ch.mix(snbam_ch.version)

		rows = snbam_ch.rows.filter{T->T.genomeId.equals(genomeId)}

		compile_ch = COMPILE([:])
		version_ch = version_ch.mix(compile_ch.version)

		cyto_ch = DOWNLOAD_CYTOBAND([:],genomeId)
                version_ch = version_ch.mix(cyto_ch.version)

		refgene_ch = DOWNLOAD_REFGENE([:],genomeId)
                version_ch = version_ch.mix(refgene_ch.version)

		prepare_ch = PREPARE_IGVREPORTS([:], genomeId, compile_ch.output, rows.map{T->T.sample+"\t"+T.bam}.collect(), vcf )
                version_ch = version_ch.mix(prepare_ch.version)

		report_ch = IGVREPORT([:], genomeId, vcf, cyto_ch.output, refgene_ch.output, prepare_ch.output.splitCsv(header:true,sep:'\t') )
                version_ch = version_ch.mix(report_ch.version)
		
		version_ch = MERGE_VERSION("IGV report", version_ch.collect())

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
	private final int MAX_CONTROLS= ${params.num_controls};
	private final int MAX_CASES= ${params.num_cases};
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
    	    final int SV_FLANKING= ${params.sv_flanking?:250};
            try(BufferedReader in=Files.newBufferedReader(bams)) {
                    in.lines().
                    	filter(T->!T.startsWith("#")).
                    	map(T->tab.split(T)).
                    	forEach(T->sample2bam.put(T[0],Paths.get(T[1])));
                    }
            final ToIntFunction<Genotype> countALT = G->(int)(G.getAlleles().stream().filter(A->!(A.isReference() || A.isNoCall())).count());
            final Comparator<Genotype> gtSorter = (A,B)->{
		if(A.isNoCall() && B.isHomRef()) return 1;
		if(A.isHomRef() && B.isNoCall()) return -1;
		if(A.isFiltered() && !B.isFiltered()) return  1;
		if(!A.isFiltered() && B.isFiltered()) return -1;
            	if(A.hasGQ() && B.hasGQ()) {
            		return Integer.compare(B.getGQ(), A.getGQ());
            		}
            	if(A.hasDP() && B.hasDP())  {
            		return Integer.compare(B.getDP(), A.getDP());
            		}
            	return A.getSampleName().compareTo(B.getSampleName());
            	};
            	
            final List<String> output_keys = Arrays.asList("title","flanking","bedpe","bams","vcf","info","cnv");
            System.out.println(String.join("\\t", output_keys));
            try(VCFIterator r = new VCFIteratorBuilder().open(System.in)) {
                long variant_id =0L;
                final VCFHeader header = r.getHeader();
                while(r.hasNext()) {
                        final VariantContext ctx = r.next();
                        final String svType = ctx.getAttributeAsString("SVTYPE","");
                        final int ctx_flanking = ctx.getAttributeAsInt("FLANKING",SV_FLANKING);
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
                            	limit(MAX_CASES).
                            	collect(Collectors.toList());
                        if(cases.isEmpty()) continue;
                        
                        final List<Genotype> controls = ctx.getGenotypes().stream().
                        	filter(G->sample2bam.containsKey(G.getSampleName())).
                        	filter(G->countALT.applyAsInt(G)==0L).
                        	sorted(gtSorter).
                        	limit(MAX_CONTROLS).
                        	collect(Collectors.toList());

                        
                        variant_id++;
                        final Set<String> sample_names = Stream.concat(cases.stream(),controls.stream()).
                        		map(G->G.getSampleName()).
                        		collect(Collectors.toSet());
                        
                        final Map<String,String> properties = new HashMap<>();
                        StringBuilder vcNameBuilder  = new StringBuilder(ctx.getContig()).append("_").append(ctx.getStart());
                        if(ctx.getEnd()!=ctx.getStart()) vcNameBuilder.append("_").append(ctx.getEnd()).append("_").append(ctx.getLengthOnReference());
                        if(!svType.isEmpty())  vcNameBuilder.append("_").append(svType);
                        final String vcName= vcNameBuilder.toString();
                        properties.put("title",vcName);
			properties.put("cnv", (!svType.isEmpty() && !svType.equals("BND"))?"true":"false");

                        
                        properties.put("bams",sample_names.stream().
                        		map(S->sample2bam.get(S).toString()).
                        		collect(Collectors.joining(" ")));
                        
                        if(svType.equals("BND")) {
                        	final List<Locatable> bnds = parseBnd(ctx,dict);
                        	if(!bnds.isEmpty()) {
                        		properties.put("flanking","100");
                        		final Path outpath = outdir.resolve(String.format("%05d",(variant_id)) + ".bedpe");
                        		properties.put("bedpe", outpath.toAbsolutePath().toString());
                        		try(PrintWriter out= new PrintWriter(Files.newBufferedWriter(outpath))) {
		                        	for(int i=0;i< bnds.size();i++) {
		                        		Locatable bnd = bnds.get(i);
		                        		out.print(ctx.getContig());
		                        		out.print("\t");
		                        		out.print(ctx.getStart()-1);
		                        		out.print("\t");
		                        		out.print(ctx.getEnd());
		                        		out.print("\t");
		                        		out.print(bnd.getContig());
		                        		out.print("\t");
		                        		out.print(bnd.getStart()-1);
		                        		out.print("\t");
		                        		out.print(bnd.getEnd());
		                        		out.print("\t");
		                        		out.println(vcName+"."+(i+1));
		                        		}
		                        	out.flush();
	                        		}
	                        	}
                        	}
                        else if(ctx.getLengthOnReference() > ctx_flanking*2) {
                        	properties.put("flanking",String.valueOf(ctx_flanking));
                        	final Path outpath = outdir.resolve(String.format("%05d",(variant_id)) + ".bedpe");
                        	properties.put("bedpe", outpath.toAbsolutePath().toString());
                        	try(PrintWriter out= new PrintWriter(Files.newBufferedWriter(outpath))) {
                    			final int chromLen = dict.getSequence(ctx.getContig()).getSequenceLength();
                        		out.print(ctx.getContig());
                        		out.print("\t");
                        		out.print(Math.max(0, ctx.getStart()-(ctx_flanking+dx5+1)));
                        		out.print("\t");
                        		out.print(Math.min(chromLen, ctx.getStart()+dx5+ctx_flanking));
                        		out.print("\t");
                        		out.print(ctx.getContig());
                        		out.print("\t");
                        		out.print(Math.max(0, ctx.getEnd()-(dx3+ctx_flanking+1)));
                        		out.print("\t");
                        		out.print(Math.min(chromLen, ctx.getEnd()+dx3+ctx_flanking));
                        		out.print("\t");
                        		out.println(vcName);
	                        	out.flush();
                        		}
                        	}
                        else {
                        	properties.put("info",ctx.getAttributes().keySet().
                        			stream().
                        			filter(K->K.equals("CSQ") || K.equals("ANN")).
                        			collect(Collectors.joining(" ")));
                        	
	                    	final VariantContextWriterBuilder vcb = new VariantContextWriterBuilder();
	                    	final Path outvcfpath = outdir.resolve(String.format("%05d",(variant_id)) + FileExtensions.COMPRESSED_VCF);
	                    	properties.put("vcf", outvcfpath.toAbsolutePath().toString());
	                    	LOG.info("write " + outvcfpath);
	                    	final VCFHeader header2 = new VCFHeader(
	                    			header.getMetaDataInInputOrder(),
	                    			header.getSampleNamesInOrder().stream().filter(S->sample_names.contains(S)).collect(Collectors.toSet())
	                    			);
	                    	
	                    	vcb.setReferenceDictionary(dict);
	                    	vcb.setIndexCreator(new TabixIndexCreator(dict, TabixFormat.VCF));
	                    	vcb.setOption(Options.INDEX_ON_THE_FLY);
	                    	
	                    	try(VariantContextWriter w = vcb. 
	                    			setOutputPath(outvcfpath).
	                    			setCreateMD5(false).
	                    			build()) {
	                    		w.writeHeader(header2);
	                    		w.add( new VariantContextBuilder(ctx).
	                    				genotypes(Stream.concat(cases.stream(), controls.stream()).collect(Collectors.toList())).
	                    				make());
	                    		}
	                        }
                        
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

process PREPARE_IGVREPORTS {
input:
	val(meta)
	val(genomeId)
	path(jar)
	val(L)
	path(vcf)
output:
	path("output.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def reference = params.genomes[genomeId].fasta
"""        
hostname 1>&2     
${moduleLoad("jvarkit bcftools")}
set -o pipefail
mkdir -p OUT

cat << EOF > jeter.list
${L.join("\n")}
EOF

bcftools view "${vcf}" |\
	java -cp ${jar}:\${JVARKIT_DIST}/vcffilterjdk.jar Minikit --reference "${reference}" --bams jeter.list --out "\${PWD}/OUT" > output.tsv


rm jeter.list
###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">prepare data for IGV reports</entry>
        <entry key="versions"></entry>
</properties>
EOF
"""
}

process IGVREPORT {
tag "${row.title}"
afterScript "rm -rf TMP"
conda "${params.conda}/IGVREPORTS"
input:
	val(meta)
	val(genomeId)
	val(vcf)
	val(cytoband)
	val(refgene)
	val(row)
output:
	path("${row.title}.html"),emit:output
	path("version.xml"),emit:version
script:
	def reference = params.genomes[genomeId].fasta
"""
hostname 1>&2
mkdir -p TMP

create_report ${row.vcf.isEmpty()?row.bedpe:row.vcf}  ${reference} \
	--ideogram "${cytoband}" \
	${row.flanking.isEmpty()?"":"--flanking ${row.flanking}"} \
	${row.info.isEmpty()?"":"--info-columns ${row.info}"} \
	--tracks ${vcf} ${row.bams} ${refgene} \
	--output TMP/${row.title}.html

mv -v "TMP/${row.title}.html" ./

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">generate IGV reports</entry>
        <entry key="versions"></entry>
</properties>
EOF
"""
}
