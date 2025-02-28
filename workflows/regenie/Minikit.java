//package com.github.lindenb.jvarkit.tools.test;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

public class Minikit extends Launcher {
	private static final String CADD_PHRED = "CADD_PHRED";
	private static final String SLIDING_WINDOW = "sliding_window";
	private static final String GNOMAD_AF = "gnomad_genome_AF_NFE";

	@Parameter(names = "-w", description = "window size")
	private int window_size = -1;
	@Parameter(names = "-a", description = "output annotation file", required = true)
	private Path annotationFileOut = null;
	@Parameter(names = "-s", description = "set list file output", required = true)
	private Path setListFileOut = null;
	@Parameter(names = "-m", description = "mask file output", required = true)
	private Path maskFileOut = null;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();

	
	private static class Variation {
		String contig;
		int pos;
		String id;
		String gene;
		String prediction;
		double score;
		double cadd;
	}

	
	private static class Prediction {
		final String name;
		final double score;
		boolean found = false;
		final Set<String> masks = new HashSet<>();

		Prediction(String name,double score) {
			this.name  = name;
			this.score  = score;
			if(name.startsWith(SLIDING_WINDOW)) {
				masks.add(this.name);
				}
			else
				{
				masks.add("ALL");
				}
			}

		void andMask(String... array) {
			for (String s : array)
				this.masks.add(s);
			}
		}

	private Prediction sliding_prediction = null;
	private final Map<String, Prediction> scores = new HashMap<>(50);

	private static class RowCodec extends AbstractDataCodec<Variation> {
		@Override
		public void encode(DataOutputStream dos, Variation variant) throws IOException {
			dos.writeUTF(variant.contig);
			dos.writeInt(variant.pos);
			dos.writeUTF(variant.id);
			dos.writeUTF(variant.gene);
			dos.writeUTF(variant.prediction);
			dos.writeDouble(variant.score);
			dos.writeDouble(variant.cadd);
		}

		@Override
		public Variation decode(DataInputStream dis) throws IOException {
			final Variation v = new Variation();
			try {
				v.contig = dis.readUTF();
			} catch (EOFException err) {
				return null;
			}
			v.pos = dis.readInt();
			v.id = dis.readUTF();
			v.gene = dis.readUTF();
			v.prediction = dis.readUTF();
			v.score = dis.readDouble();
			v.cadd = dis.readDouble();
			return v;
			}

		@Override
		public RowCodec clone() {
			return new RowCodec();
		}
	}

	private Prediction makeScore(String pred, double score) {
		final Prediction p = new Prediction(pred,score);
		this.scores.put(pred, p);
		return p;
	}

	public Minikit() {

	}

	private String fixContig(final String ctg) {
		if(ctg.equals("X") || ctg.equals("chrX")) return "23";
		if(ctg.equals("Y") || ctg.equals("chrY")) return "24";
		if(ctg.startsWith("chr")) return ctg.substring(3);
		return ctg;
		}

	private void dump(final SortingCollection<Variation> sorter, final AnnPredictionParser annParser, final VariantContext ctx) throws Exception {
		if(ctx.hasAttribute(GNOMAD_AF) && ctx.getAttributeAsDouble(GNOMAD_AF, 0.0) >= 0.1) return;
		if(ctx.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY) && ctx.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, 0.0) >= 0.1) return;
		
		final int  win_pos =  this.window_size<1 ? -1:
				1 + (((int)(ctx.getStart()/(double)this.window_size)) * this.window_size)
				;
		Double cadd_phred = null;
		if (ctx.hasAttribute(CADD_PHRED)) {
				String s = ctx.getAttributeAsString(CADD_PHRED, ".");
				if (!(s.equals(".") || StringUtils.isBlank(s))) {
					cadd_phred = Double.valueOf(s);
				}
			}
		
		Double best_score = null;
		Prediction best_pred = null;
		String best_gene = null;
		
		final String altstr = ctx.getAlternateAllele(0).getDisplayString();
		final List<AnnPredictionParser.AnnPrediction> predictions = annParser.getPredictions(ctx);
			
		for(AnnPredictionParser.AnnPrediction pred:predictions) {
			String gene_name = pred.getGeneName();
			if(gene_name.isEmpty()|| gene_name.equals(".")) {
				gene_name= pred.getGeneId();
				}
			if(gene_name.isEmpty() || gene_name.equals(".")) continue;
			if(!pred.getAllele().equalsIgnoreCase(altstr)) continue;
			for(String pred_key : pred.getSOTermsStrings()) {
				if(pred_key.equals("intergenic_region")) continue;
				final Prediction p = scores.getOrDefault(pred_key, null);
				if (p == null) throw new IOException("undefined prediction key "+pred_key);
				if (best_score == null || best_score.compareTo(p.score) < 0) {
					best_score = p.score;
					best_pred = p;
					best_gene = gene_name;
					}
				}
			}
					

		if (best_score != null) {
			best_pred.found=true;
			
			final Variation v = new Variation();
			v.contig = fixContig(ctx.getContig());
			v.pos = ctx.getStart();
			v.id = ctx.getID();
			v.gene = best_gene;
			v.prediction = best_pred.name;
			v.score = best_score;
			v.cadd = (cadd_phred==null?0.0:cadd_phred.doubleValue());
			sorter.add(v);
			}

		if(win_pos>0) 
			{
			sliding_prediction.found = true;
			
			final Variation v = new Variation();
			v.contig = fixContig(ctx.getContig());
			v.pos = ctx.getStart();
			v.id = ctx.getID();
			v.gene = ctx.getContig()+ "_" + (win_pos) + "_" + (win_pos - 1 + this.window_size) ;
			v.prediction = sliding_prediction.name;
			v.score = sliding_prediction.score;
			v.cadd = (cadd_phred==null?0.0:cadd_phred.doubleValue());
			sorter.add(v);
			}
		}

	private static int compareX(Variation t1, Variation t2, int level) {
		int i = t1.contig.compareTo(t2.contig);// CHROMOSOME
		if (i != 0)
			return i;
		i = t1.gene.compareTo(t2.gene);// GENE
		if (i != 0)
			return i;
		if (level == 1)
			return 0;
		return Integer.compare(t1.pos, t2.pos);
	}

	@Override
	public int doWork(List<String> args) {
		try {
			if(this.window_size>0) {
				this.sliding_prediction = makeScore(SLIDING_WINDOW+"_"+this.window_size , 0.1);
				}
			else
				{
				this.sliding_prediction = null;
				makeScore("3_prime_UTR_variant", 0.1).andMask("UTR", "UTR3");
				makeScore("5_prime_UTR_premature_start_codon_gain_variant", 0.2).andMask("UTR", "UTR5");
				makeScore("5_prime_UTR_truncation", 0.5).andMask("UTR", "UTR5");
				makeScore("3_prime_UTR_truncation", 0.5).andMask("UTR", "UTR3");
				makeScore("5_prime_UTR_variant", 0.2).andMask("UTR", "UTR5");
				makeScore("bidirectional_gene_fusion", 1.0);
				makeScore("conservative_inframe_deletion", 0.1).andMask("protein_altering");
				makeScore("conservative_inframe_insertion", 0.3).andMask("protein_altering");;
				makeScore("disruptive_inframe_deletion", 0.2).andMask("protein_altering");
				makeScore("disruptive_inframe_insertion", 0.2).andMask("protein_altering");
				makeScore("downstream_gene_variant", 0.1).andMask("downstream", "updownstream");
				makeScore("exon_loss_variant", 1.0).andMask("protein_altering");
				makeScore("frameshift_variant", 0.4);
				makeScore("exon_loss_variant", 1.0).andMask("protein_altering");
				makeScore("gene_fusion", 0.9).andMask("protein_altering");
				makeScore("intergenic_region", 0.001);
				makeScore("intragenic_variant", 0.01);
				makeScore("initiator_codon_variant",0.5).andMask("protein_altering");
				makeScore("intron_variant", 0.05).andMask("intronic");
				makeScore("missense_variant", 0.9).andMask("protein_altering");
				makeScore("non_coding_transcript_exon_variant", 0.1).andMask("non_coding");
				makeScore("non_coding_transcript_variant", 0.1).andMask("non_coding");
				makeScore("splice_acceptor_variant", 0.5).andMask("protein_altering", "splice");
				makeScore("splice_donor_variant", 0.5).andMask("protein_altering", "splice");
				makeScore("splice_region_variant", 0.5).andMask("protein_altering", "splice");
				makeScore("start_retained_variant",0.1).andMask("synonymous");
				makeScore("start_lost", 0.6).andMask("protein_altering");
				makeScore("stop_gained", 0.9).andMask("protein_altering");
				makeScore("stop_lost", 0.6).andMask("protein_altering");
				makeScore("stop_retained_variant", 0.2).andMask("synonymous");
				makeScore("synonymous_variant", 0.1).andMask("synonymous");
				makeScore("upstream_gene_variant", 0.1).andMask("upstream", "updownstream");
				}
			final SortingCollection<Variation> sorter = SortingCollection.newInstance(Variation.class, new RowCodec(),
					(A, B) -> compareX(A, B, 2), writingSortingCollection.getMaxRecordsInRam(),
					writingSortingCollection.getTmpPaths());
			sorter.setDestructiveIteration(true);


			try (VCFIterator iter = new VCFIteratorBuilder().open(System.in)) {
				final VCFHeader vcfHeader = iter.getHeader();
				final AnnPredictionParser annParser = new AnnPredictionParserFactory(vcfHeader).get();
				if(!annParser.isValid()) throw new IOException("cannot create ANN parser");

				while (iter.hasNext()) {
					final VariantContext vc = iter.next();
					if (vc.getNAlleles() != 2)
						throw new IOException(vc.getContig() + ":" + vc.getStart() + ":" + vc.getAlleles());
					dump(sorter, annParser, vc);
				} // end while
			}
			sorter.doneAdding();

			try (CloseableIterator<Variation> iter0 = sorter.iterator()) {
				try (EqualRangeIterator<Variation> iter = new EqualRangeIterator<>(iter0, (A, B) -> compareX(A, B, 1))) {

					try (PrintWriter annotOut = IOUtils.openPathForPrintWriter(this.annotationFileOut)) {
						try (PrintWriter setFileOut = IOUtils.openPathForPrintWriter(this.setListFileOut)) {

							while (iter.hasNext()) {
								final List<Variation> gene_variants = iter.next();
								final Variation first = gene_variants.get(0);

								for (Variation v : gene_variants) {
									annotOut.print(v.id);
									annotOut.print(" ");
									annotOut.print(v.gene);
									annotOut.print(" ");
									annotOut.print(v.prediction);
									annotOut.print(" ");
									annotOut.print(v.score);
									annotOut.print(" ");
									annotOut.print(v.cadd);
									annotOut.println();
									}
								setFileOut.print(first.gene);// gene
								setFileOut.print("\t");
								setFileOut.print(first.contig);// contig
								setFileOut.print("\t");
								setFileOut.print((int) gene_variants.stream().mapToInt(it -> it.pos).average().getAsDouble());
								setFileOut.print("\t");
								setFileOut.print(gene_variants.stream().map(it -> it.id).collect(Collectors.joining(",")));
								setFileOut.println();

							}
							setFileOut.flush();
						}
						annotOut.flush();
					}
				}

			}

			sorter.cleanup();
			
			try (PrintWriter maskOut = IOUtils.openPathForPrintWriter(this.maskFileOut)) {
				final Set<String> seen = scores.values().stream().filter(P -> P.found).flatMap(S -> S.masks.stream())
						.collect(Collectors.toSet());
				for (final String mask_name : seen) {
					maskOut.print(mask_name);
					maskOut.print("\t");
					maskOut.print(scores.values().stream().filter(P -> P.masks.contains(mask_name)).
							map(P -> P.name).collect(Collectors.joining(",")));
					maskOut.println();
				}
			}
			return 0;
		} catch (Throwable err) {
			err.printStackTrace();
			return -1;
		}
	}

	public static void main(final String[] args) {
		new Minikit().instanceMainWithExit(args);
	}
}
