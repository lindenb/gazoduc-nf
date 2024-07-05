/*
The MIT License (MIT)

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
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.function.ToDoubleFunction;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.gatk.GATKConstants;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.iterator.EqualIterator;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.stats.FisherCasesControls;
import com.github.lindenb.jvarkit.pedigree.CasesControls;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.predictions.GeneExtractorFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.GeneExtractorFactory.GeneExtractor;

import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;
import htsjdk.variant.vcf.VCFRecordCodec;


public class Minikit extends Launcher {
	private static final Logger LOG = Logger.build( Minikit.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();
	@ParametersDelegate
	private CasesControls casesControls = new CasesControls();
	@Parameter(names={"--header"},description = "insert this file in the output header")
	private Path headerFile = null;
	@Parameter(names={"--body"},description = "insert this file after each output")
	private Path bodyFile = null;	
	@Parameter(names={"--minGQsingleton"},description = "min singleton GQ")
	private int minGQsingleton=90;
	@Parameter(names={"--minRatioSingleton"},description = "min singleton ration AD")
	private double minRatioSingleton=0.3;
	@Parameter(names={"--"+GATKConstants.QD_KEY},description = "GATK "+GATKConstants.QD_KEY)
	private double gatk_QD=-1.0;
	@Parameter(names={"--"+GATKConstants.FS_KEY},description = "GATK "+GATKConstants.FS_KEY)
	private double gatk_FS=-1.0;
	@Parameter(names={"--"+GATKConstants.SOR_KEY},description = "GATK "+GATKConstants.SOR_KEY)
	private double gatk_SOR=-1.0;
	@Parameter(names={"--"+GATKConstants.MQ_KEY},description = "GATK "+GATKConstants.MQ_KEY)
	private double gatk_MQ=1000;
	@Parameter(names={"--"+GATKConstants.MQRankSum_KEY},description = "GATK "+GATKConstants.MQRankSum_KEY)
	private double gatk_MQRankSum=1000;
	@Parameter(names={"--"+GATKConstants.ReadPosRankSum_KEY},description = "GATK "+GATKConstants.ReadPosRankSum_KEY)
	private double gatk_ReadPosRankSum=1000;
	@Parameter(names={"--minDP"},description = "minDP")
	private int minDP=0;
	@Parameter(names={"--maxDP"},description = "maxDP")
	private int maxDP=300;
	@Parameter(names={"--lowGQ"},description = "lowGQ")
	private int lowGQ=60;
	@Parameter(names={"--f_missing"},description = "fraction missing")
	private double f_missing=0.05;
	@Parameter(names={"--gtf"},description = "gtf indexed file")
	private String gtfPath=null;


	private TabixReader gtfReader = null;
	
	
	private static class Key implements Comparable<Key>{
		String contig; // eg. chr1
		String splitter; // eg. SNpeff/transcript
		String key;// e.g ENST00000001
		String geneName; // e.g. SCN5A
		
		@Override
		public int hashCode() {
			int i=contig.hashCode();
			i=i%31*splitter.hashCode();
			i=i%31*key.hashCode();
			return i;
			}
		
		@Override
		public int compareTo(Key o) {
			int i= this.contig.compareTo(o.contig);
			if(i!=0) return i;
			i= this.splitter.compareTo(o.splitter);
			if(i!=0) return i;
			i= this.key.compareTo(o.key);
			if(i!=0) return i;
			// do not use geneName
			return 0;
			}
		@Override
		public boolean equals(Object obj) {
			if(obj==this) return true;
			if(obj==null || !(obj instanceof Key)) return false;
			return compareTo(Key.class.cast(obj))==0;
			}
		@Override
		public String toString() {
			return contig+" "+splitter+" "+key+" "+geneName;
			}
		}
	
private interface Splitter extends Function<VariantContext,List<Key>> {
}

private class GeneSplitter implements Splitter {
	final List<GeneExtractor> extractors ;
	public GeneSplitter(VCFHeader header) {
		GeneExtractorFactory geneExtractorFactory= new GeneExtractorFactory(header);
		this.extractors = geneExtractorFactory.getAllExtractors();
		}
	@Override
	public List<Key> apply(VariantContext vc) {
		final Set<Key> keys = new HashSet<>();
		for(GeneExtractor ex: this.extractors) {
			for(GeneExtractorFactory.KeyAndGene kgx:ex.apply(vc).keySet()) {
				final Key k= new  Key();
				k.contig=vc.getContig();
				k.splitter = kgx.getMethod();
				k.key = kgx.getKey();
				k.geneName = StringUtils.ifBlank(kgx.getGene(),".");
				keys.add(k);
				}
			}
		return new ArrayList<>(keys);
		}
}


private static class KeyAndVariant {
	Key key;
	VariantContext ctx;
	KeyAndVariant(final Key key, final VariantContext ctx) {
		this.key = key;
		this.ctx = ctx;
		}

	public int compareTo2(KeyAndVariant o) {
		int i= this.key.compareTo(o.key);
		if(i!=0) return i;
		i = this.ctx.getContig().compareTo(o.ctx.getContig());
		if(i!=0) return i;
		i = Integer.compare(this.ctx.getStart(),o.ctx.getStart());
		if(i!=0) return i;
		i = this.ctx.getReference().compareTo(o.ctx.getReference());
		if(i!=0) return i;
		int x=1;
		while(x<this.ctx.getAlleles().size() && x<o.ctx.getAlleles().size()) {
			i = this.ctx.getAlleles().get(x).compareTo(o.ctx.getAlleles().get(x));
			if(i!=0) return i;
			++x;
			}
		return 0;
		}
	}

private static class KeyAndVariantCodec implements SortingCollection.Codec<KeyAndVariant> {
	final VCFHeader header;
	private final VCFRecordCodec vcfCodec;
	private final BinaryCodec binaryCodecR = new BinaryCodec();
	private final BinaryCodec binaryCodecW = new BinaryCodec();
 
	KeyAndVariantCodec(VCFHeader header) {
		this.header=header;
		this.vcfCodec = new VCFRecordCodec(header);
		}
	@Override
	public void setInputStream(InputStream is) {
		this.vcfCodec.setInputStream(is);
		this.binaryCodecR.setInputStream(is);
		}
	@Override
	public void setOutputStream(OutputStream os) {
		this.vcfCodec.setOutputStream(os);
		this.binaryCodecW.setOutputStream(os);
		}
	@Override
	public KeyAndVariantCodec clone() {
		return new KeyAndVariantCodec(this.header);
		}
	@Override
	public KeyAndVariant decode() {
		final Key key = new Key();
		try {
			key.contig = this.binaryCodecR.readNullTerminatedString();
			}
		catch(Throwable err) {
			return null;
			}
		key.splitter = this.binaryCodecR.readNullTerminatedString();
		key.key = this.binaryCodecR.readNullTerminatedString();
		key.geneName = this.binaryCodecR.readNullTerminatedString();
		VariantContext ctx = this.vcfCodec.decode();
		return new KeyAndVariant(key, ctx);
		}
	private void writeNullTerminatedString(String s) {
		this.binaryCodecW.writeString(s, false, true);
		}
	@Override
	public void encode(KeyAndVariant o) {
		writeNullTerminatedString(o.key.contig);
		writeNullTerminatedString(o.key.splitter);
		writeNullTerminatedString(o.key.key);
		writeNullTerminatedString(o.key.geneName);
		this.vcfCodec.encode(o.ctx);
		}
	}

private boolean accept(final VariantContext variant) {
	if(variant.isFiltered()) {
		return false;
	}
	

	final int count_alt = (int)variant.getGenotypes().stream().
        filter(g->g.hasAltAllele()).
        count();

	if(count_alt==0) return false;


	
	
	if(this.f_missing  >= 0 ) {
	        /** missing */
	        final double  n_missing = variant.getGenotypes().stream().
	                filter(G->G.isNoCall() || (G.hasDP() && G.getDP()==0)).
	                count();

	        if(n_missing/variant.getNSamples() > this.f_missing) return false;
	        }


	
	/** low DP */
	final double dp= variant.getGenotypes().stream().
	        filter(G->G.isCalled() && G.hasDP()).
	        mapToInt(G->G.getDP()).average().orElse(this.minDP);

	if(dp < this.minDP  || dp > this.maxDP) return false;

	
	final int count_alt_filtered=(int)variant.getGenotypes().stream().
			filter(G->G.hasAltAllele() && G.isFiltered()).
			count();
	if(count_alt_filtered/count_alt >= 0.1) {
		return false;
		}
	
	

	/** low GQ */
	final double count_low_gq = variant.getGenotypes().stream().
	        filter(G->G.hasAltAllele() && G.hasGQ()).
	        filter(g->g.getGQ()< this.lowGQ).
	        count();
	
	if(count_low_gq/count_alt >= 0.25) {
	        return false;
	        }

	
	
	if(variant.hasAttribute(GATKConstants.QD_KEY) &&
		variant.getAttributeAsDouble(GATKConstants.QD_KEY,1000) < this.gatk_QD) {
		return false;
		}
	
	if(variant.hasAttribute(GATKConstants.FS_KEY) && 
			variant.getAttributeAsDouble(GATKConstants.FS_KEY,0) > this.gatk_FS) {
		return false;
		}
	
	if(variant.hasAttribute(GATKConstants.SOR_KEY) && 
			variant.getAttributeAsDouble(GATKConstants.SOR_KEY,0) > this.gatk_SOR) {
		return false;
		}
	
	if(variant.hasAttribute(GATKConstants.MQ_KEY) && 
			variant.getAttributeAsDouble(GATKConstants.MQ_KEY,1000) < this.gatk_MQ) {
		return false;
		}
	
	if(variant.hasAttribute(GATKConstants.MQRankSum_KEY) && 
			variant.getAttributeAsDouble(GATKConstants.MQRankSum_KEY,1000) < this.gatk_MQRankSum) {
		return false;
		}
	
	if(variant.hasAttribute(GATKConstants.ReadPosRankSum_KEY) && 
			variant.getAttributeAsDouble(GATKConstants.ReadPosRankSum_KEY,1000) < this.gatk_ReadPosRankSum) {
		return false;
		}
	

	
	Genotype singleton=null;
	for(final Genotype g: variant.getGenotypes()) {
	        if(g.hasAltAllele()) {
                if(singleton!=null) {
                	singleton=null;
                	break;
                	}
                singleton=g;
                }
	        }
	if(singleton!=null) {
		if(singleton.isFiltered()) return true;
		if(singleton.isHet() && singleton.hasGQ() && singleton.getGQ()< minGQsingleton) {
	        return false;
	        }
		if(singleton.hasAD() && singleton.isHet() && singleton.getAD().length==2)
	        {
			final  int array[]=singleton.getAD();
	        final double r= array[1]/(double)(array[0]+array[1]);
	        if(r< this.minRatioSingleton || r>(1.0 - minRatioSingleton)) return false;
	        }
		}
	return true;
	}

@Override
public int doWork(List<String> args) {
	try {
		if(!StringUtils.isBlank(this.gtfPath)) {
			this.gtfReader = new TabixReader(this.gtfPath);
			}
		
		try(VCFIterator iter = new VCFIteratorBuilder().open(System.in)) {
			final VCFHeader header = iter.getHeader();

			this.casesControls.load();
			this.casesControls.retain(header);
			this.casesControls.checkHaveCasesControls();
			final Set<String> all_samples = casesControls.getAll();

            final UnaryOperator<String> noChr = S->S.startsWith("chr")?S.substring(3):S;

			final ToDoubleFunction<VariantContext> toAF=VC->{
				double ac=0;
				double an=0;
				for(String sn:all_samples) {
					final Genotype gt = VC.getGenotype(sn);
					if(gt==null) continue;
					ac+= (int)gt.getAlleles().stream().filter(A->!(A.isReference() || A.isNoCall())).count();
					an+= gt.getPloidy();
					}
				if(an==0.0) return 0.0;
				return ac/an;
				};
			
			final SortingCollection<KeyAndVariant> sorter= SortingCollection.newInstance(
					KeyAndVariant.class,
					new KeyAndVariantCodec(header),
					(A,B)->A.compareTo2(B),
					writingSortingCollection.getMaxRecordsInRam(),
					writingSortingCollection.getTmpPaths()
					);
			
			final Splitter splitter= new GeneSplitter(header);
			while(iter.hasNext()) {
				final VariantContext ctx=iter.next();
				if(!accept(ctx)) continue;
				for(Key key: splitter.apply(ctx)) {
					sorter.add(new KeyAndVariant(key,ctx));
					}
				}
			sorter.doneAdding();
			try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				if(headerFile!=null) {
					IOUtils.copyTo(headerFile, pw);
					pw.println();
					}
				pw.println("# samples ( 0: unaffected 1:affected)");

				pw.print("population <- data.frame(family=c(");
				pw.print(all_samples.stream().map(S->StringUtils.doubleQuote(S)).collect(Collectors.joining(",")));
				pw.print("),name=c(");
				pw.print(all_samples.stream().map(S->StringUtils.doubleQuote(S)).collect(Collectors.joining(",")));
				pw.print("),status=c(");
				pw.print(all_samples.stream().map(S->casesControls.isControl(S)?"0":"1").collect(Collectors.joining(",")));
				pw.println("))");

				try(CloseableIterator<KeyAndVariant> iter0= sorter.iterator()) {
					final EqualIterator<KeyAndVariant> iter1 = new EqualIterator<>(iter0,(A,B)->A.key.compareTo(B.key));
					while(iter1.hasNext()) {
						final List<KeyAndVariant> array = iter1.next();
						final Key key = array.get(0).key;
						final List<VariantContext> variants = array.stream().map(it->it.ctx).collect(Collectors.toList());
						
						// no Case have an ALT, skip
						if(variants.stream().flatMap(V->V.getGenotypes().stream()).
							filter(G->G.hasAltAllele()).
							noneMatch(G->casesControls.isCase(G))) continue;
						
						final FisherCasesControls fisherCasesControls = new FisherCasesControls(casesControls);

						pw.println("contig <- "+ StringUtils.doubleQuote(key.contig));
						pw.println("splitter <- "+ StringUtils.doubleQuote(key.splitter));
						pw.println("gene.key <- "+ StringUtils.doubleQuote(key.key));
						pw.println("gene.name <- "+ StringUtils.doubleQuote(key.geneName));
						pw.println("n.variants <- "+ variants.size());
						pw.println("gene.start <- "+ variants.stream().mapToInt(V->V.getStart()).min().orElse(-1));
						pw.println("gene.end <- "+ variants.stream().mapToInt(V->V.getStart()).max().orElse(-1));
						pw.println("variants.str <- "+ StringUtils.doubleQuote(variants.stream().map(V->noChr.apply(V.getContig())+":"+V.getStart()).collect(Collectors.joining(";"))));

						pw.print("genotypes <- c(");
							for(int i=0;i< variants.size();i++) {
							final VariantContext ctx =  variants.get(i);
							if(i>0) pw.println(",");
							pw.print(all_samples.stream().
								map(S->ctx.getGenotype(S)).
								map(genotype->{
									fisherCasesControls.accept(genotype);
									final int n= (int)genotype.getAlleles().stream().filter(A->!(A.isReference() || A.isNoCall())).count();
									switch(n) {
										case 0: return "0";
										case 1: return "1";
										default: return "2";
										}
									}).
								collect(Collectors.joining(","))
								);
							}// end reading vcf
						pw.println(")");
						pw.flush();

						pw.println("fisher.pvalue <- "+ fisherCasesControls.getAsDouble());
						pw.println("CASES_ALT_COUNT <- "+ fisherCasesControls.getCasesAltCount());
						pw.println("CASES_REF_COUNT <- "+ fisherCasesControls.getCasesRefCount());
						pw.println("CTRLS_ALT_COUNT <- "+ fisherCasesControls.getControlsAltCount());
						pw.println("CTRLS_REF_COUNT <- "+ fisherCasesControls.getControlsRefCount());

						
						pw.print("# variants.");
						pw.println();
						pw.print("variants <- data.frame(chrom=c(");
						pw.print(variants.stream().map(vc->StringUtils.doubleQuote(vc.getContig())).collect(Collectors.joining(",")));
						pw.print("),chromStart=c(");
						pw.print(variants.stream().map(vc->String.valueOf(vc.getStart())).collect(Collectors.joining(",")));
						pw.print("),chromEnd=c(");
						pw.print(variants.stream().map(vc->String.valueOf(vc.getEnd())).collect(Collectors.joining(",")));
						pw.print("),refAllele=c(");
						pw.print(variants.stream().map(vc->StringUtils.doubleQuote(vc.getReference().getDisplayString())).collect(Collectors.joining(",")));
						pw.print("),altAllele=c(");
						pw.print(variants.stream().map( vc->StringUtils.doubleQuote(vc.getAlternateAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining(",")))).collect(Collectors.joining(",")));
						pw.print("),AF=c(");
						pw.print(variants.stream().mapToDouble(toAF).mapToObj(AF->String.valueOf(AF)).collect(Collectors.joining(",")));
						pw.print(")");
						pw.println(")");

						if(bodyFile!=null) {
							IOUtils.copyTo(bodyFile, pw);
							pw.println();
							}
						}
					iter1.close();
					}
				pw.flush();
				}
			}
		return 0;
		}
	catch(final Throwable err) {
		LOG.error(err);
		return -1;
		}
	finally {
		if(gtfReader!=null) gtfReader.close();
		}
	}

public static void main(final String[] args) {
	new Minikit().instanceMain(args);
	}
}
