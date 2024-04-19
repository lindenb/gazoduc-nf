import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.MergingIterator;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public class Minikit {
		private static class Variant implements Comparable<Variant> {
		int pos;
		List<Allele> alleles;
		String id;
		int priority=0;
		@Override
		public int compareTo(final Variant o) {
			int i = Integer.compare(this.pos, o.pos);
			if(i!=0) return i;
			i = this.alleles.get(0).compareTo(o.alleles.get(0));
			return i;
			}
		@Override
		public boolean equals(Object o) {
			if(o==this) return true;
			if(o==null || !(o instanceof Variant)) return false;
			return this.compareTo(Variant.class.cast(o))==0;
			}
		@Override
		public int hashCode() {
			return Integer.hashCode(pos);
			}

		@Override
		public String toString() {
			return String.valueOf(pos)+":"+id+":"+this.alleles.get(0);
		}
	}
	
	private static class VCFSource {
		final int priority;
		final Path filePath;
		VCFSource(final int priority,final Path path){
			this.priority = priority;
			this.filePath = path;
			}
	    }
	    
	private static class VcfScanner extends AbstractCloseableIterator<Variant> {
        final VCFSource source;
		final VCFFileReader reader;
		final PeekableIterator<VariantContext> iter;
		final List<Variant> stack  = new ArrayList<>();

		VcfScanner(VCFSource src,final SAMSequenceRecord ssr)   {
            this.source= src;
		    this.reader = new VCFFileReader(src.filePath, true);
		    this.iter = new PeekableIterator<>(this.reader.query(ssr));
			}
		
		private Variant convert(final VariantContext ctx) {
	 		final Variant variant = new Variant();
            variant.id = ctx.getID();
            variant.alleles = ctx.getAlleles();
            variant.pos = ctx.getStart();
            variant.priority = source.priority;
			return variant;
			}

		@Override
		protected Variant advance() {
			if(!stack.isEmpty()) {
				return stack.remove(0);
				}
			if(iter.hasNext()) {
				final VariantContext ctx = iter.next();
				final Variant variant = convert(ctx);
				stack.add(variant);
				while(iter.hasNext()) {
					final VariantContext v2 = iter.peek();
					if(v2.getStart()!=variant.pos) {
						break;
						}
					final Variant variant2 = convert(iter.next());//consumme
					stack.add(variant2);
					}
				}
			if(!stack.isEmpty()) {
				if(stack.size()>1) Collections.sort(this.stack);
				return stack.remove(0);
				}
			return null;
			}
		
		@Override
		public void close()  {
			iter.close();
			try { this.reader.close();} catch(Throwable err) {}
			}
	}


	
	public int doWork(final String[] args) {
		try {
			final Path fasta = Paths.get(args[0]);
			final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(fasta);
			
			final List<VCFSource> sources = Files.newBufferedReader(Paths.get(args[1])).lines().
				map(L->{
					final String[] tokens = L.split("[\t]");
					return new VCFSource(Integer.parseInt(tokens[0]),Paths.get(tokens[1]));
					}).collect(Collectors.toList());
			
			final VariantContextWriterBuilder vcwb = new VariantContextWriterBuilder();
			vcwb.clearOptions();
			vcwb.setReferenceDictionary(dict);
			vcwb.setOutputVCFStream(System.out);
			
			final Pattern rsIdPattern = Pattern.compile("[rR][sS][0-9]+");
			try(VariantContextWriter w=vcwb.build()) {
				final VCFHeader header = new VCFHeader();
				header.setSequenceDictionary(dict);
				w.writeHeader(header);
				for(final SAMSequenceRecord ssr:dict.getSequences()) {
					final List<CloseableIterator<Variant>> iterators  = sources.stream().map(SRC->new VcfScanner(SRC,ssr)).collect(Collectors.toList());
					final MergingIterator<Variant> merger = new MergingIterator<>((A, B)->Integer.compare(A.pos,B.pos),iterators);
					final EqualRangeIterator<Variant> equal_range = new EqualRangeIterator<>(merger, (A, B)->Integer.compare(A.pos,B.pos));
					
					while(equal_range.hasNext()) {
						final List<Variant> variants0 = equal_range.next();
						//all references in this set of variant
						final Set<Allele> all_refs = variants0.stream().map(V->V.alleles.get(0)).collect(Collectors.toCollection(TreeSet::new));
						// loop over each ref
						for(final Allele ref_allele : all_refs) {
							// all variants with this ref allele
							final List<Variant> variants = variants0.stream().filter(V->V.alleles.get(0).equals(ref_allele)).collect(Collectors.toList());
							final Variant variant = variants.get(0);
							final Set<Allele> altSet = new HashSet<>();
							variants.stream().forEach(V->altSet.addAll(V.alleles.subList(1, V.alleles.size())));
							altSet.remove(Allele.SPAN_DEL);
							final List<Allele> alleles = new ArrayList<>(1+altSet.size());
							alleles.add(ref_allele);
							alleles.addAll(altSet);
							
							
							final  VariantContextBuilder vcb = new VariantContextBuilder(null, ssr.getContig(), variant.pos,variant.pos+variant.alleles.get(0).length()-1, alleles);
							String id = variants.
									stream().
									filter(F->Arrays.stream(CharSplitter.SEMICOLON.split(F.id)).// vcf spec:  Semicolon-separated list of unique identifiers where available
											anyMatch(S->rsIdPattern.matcher(S).matches())).
									map(F->F.id).
									findFirst().
									orElse(null);
							// not an rs ID
							if(StringUtils.isBlank(id)) id = variants.stream().
									sorted((A,B)->Integer.compare(A.priority, B.priority)).
									findFirst().
									map(F->F.id).
									orElse(null);
							vcb.id(variant.id);
							w.add(vcb.make());
							}
						}
				
					equal_range.close();
					merger.close();
        			iterators.stream().forEach(SRC->SRC.close());
					}
				}
			return 0;
		} catch(final Throwable err) {
			err.printStackTrace();
			return -1;
		}
	}
	
	public static void main(final String[] args) {
		new Minikit().doWork(args);
	}

}

