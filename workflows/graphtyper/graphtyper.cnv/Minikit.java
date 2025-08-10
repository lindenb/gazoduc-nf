import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

public class Minikit extends Launcher {
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputDir = null;
	@Parameter(names={"-L","--length"},description="max-length")
	private int max_len = 100_000;
	
	private int ID_GENERATOR=0;
	
private class Batch {
	final List<VariantContext> variants =new ArrayList<>();

	long getLengthOnReference() {
		return variants.stream().mapToLong(V->V.getLengthOnReference()).sum();
		}
	void save(VCFHeader header) throws IOException {
		final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
		final File filename = new File(outputDir,String.format("split.%d.L%s.N%d.vcf.gz", 
				(++ID_GENERATOR),
				String.valueOf(getLengthOnReference()),
				variants.size()
				));
		final TabixIndexCreator tbi=new TabixIndexCreator(dict, TabixFormat.VCF);
		final VariantContextWriterBuilder vcw= new VariantContextWriterBuilder();
		vcw.setIndexCreator(tbi);
		vcw.setOutputFile(filename);
		vcw.setOutputFileType(VariantContextWriterBuilder.determineOutputTypeFromFile(filename));
		vcw.setReferenceDictionary(dict);
		vcw.setOption(Options.INDEX_ON_THE_FLY);
		vcw.unsetBuffering();
		try( VariantContextWriter  w=vcw.build()) {
			w.writeHeader(header);
			for(VariantContext vc: this.variants) {
				w.add(vc);
				}
			}
		}
	}

@Override
public int doWork(List<String> args) {

	try {
		final Random rand = new Random(0L);
		try(VCFIterator iter = new VCFIteratorBuilder().open(System.in)) {
			final VCFHeader header = iter.getHeader();
			final List<Batch> batches = new ArrayList<>();
			for(;;) {
				final VariantContext ctx=iter.hasNext()?iter.next():null;
				
				if(ctx==null) {
					for(Batch b:batches) {
						b.save(header);
						}
					break;
					}
				
				final int len = ctx.getLengthOnReference();

				if(len >= max_len) {
					final Batch b=new Batch();
					b.variants.add(ctx);
					b.save(header);
					continue;
					}
				
			
				final List<Batch> copy=new ArrayList<>(batches);
				copy.removeIf(B->B.getLengthOnReference()+len > max_len);
				if(!copy.isEmpty()) {
					Collections.shuffle(copy,rand);
					copy.get(0).variants.add(ctx);
					}
				else
					{
					final Batch b=new Batch();
					b.variants.add(ctx);
					batches.add(b);
					}
				int i=0;
				while(i< batches.size()) {
					if( (double)batches.get(i).getLengthOnReference() >= max_len*0.9 || batches.get(i).variants.size()>100) {
						batches.get(i).save(header);
						batches.remove(i);
						}
					else
						{
						i++;
						}
					}
				}
			}
	   
		return 0;
		}
	catch(final Throwable err) {
		err.printStackTrace();
		return -1;
		}
	}

public static void main(final String[] args) {
	new Minikit().instanceMain(args);
	}
}
