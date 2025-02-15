//package com.github.lindenb.jvarkit.tools.test;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

public class Minikit extends Launcher {
	private final String ANN = "ANN";
	private final String CADD_PHRED = "CADD_PHRED";

	@Parameter(names = "-a", description = "output annotation file", required = true)
	private Path annotationFileOut = null;
	@Parameter(names = "-s", description = "set list file output", required = true)
	private Path setListFileOut = null;
	@Parameter(names = "-m", description = "mask file output", required = true)
	private Path maskFileOut = null;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();

	private static class Prediction {
		String name;
		double score;
		boolean found = false;
		Set<String> masks = new HashSet<>();

		Prediction() {
			masks.add("ALL");
		}

		void andMask(String... array) {
			for (String s : array)
				this.masks.add(s);
		}
	}

	private final Map<String, Prediction> scores = new HashMap<>(50);

	private static class RowCodec extends AbstractDataCodec<String> {
		@Override
		public void encode(DataOutputStream dos, String line) throws IOException {
			dos.writeUTF(line);
		}

		@Override
		public String decode(DataInputStream dis) throws IOException {
			try {
				return dis.readUTF();
			} catch (EOFException err) {
				return null;
			}
		}

		@Override
		public RowCodec clone() {
			return new RowCodec();
		}
	}

	private Prediction makeScore(String pred, double score) {
		final Prediction p = new Prediction();
		p.name = pred;
		p.score = score;
		this.scores.put(pred, p);
		return p;
	}

	public Minikit() {
		makeScore("3_prime_UTR_variant", 0.1).andMask("UTR", "UTR3");
		makeScore("5_prime_UTR_premature_start_codon_gain_variant", 0.2).andMask("UTR", "UTR5");
		makeScore("5_prime_UTR_truncation", 0.5).andMask("UTR", "UTR5");
		;
		makeScore("5_prime_UTR_variant", 0.2).andMask("UTR", "UTR5");
		;
		makeScore("bidirectional_gene_fusion", 1.0);
		makeScore("conservative_inframe_deletion", 0.1);
		makeScore("conservative_inframe_insertion", 0.3);
		makeScore("disruptive_inframe_deletion", 0.2).andMask("protein_altering");
		makeScore("disruptive_inframe_insertion", 0.2).andMask("protein_altering");
		makeScore("downstream_gene_variant", 0.1).andMask("downstream", "updownstream");
		makeScore("exon_loss_variant", 1.0).andMask("protein_altering");
		makeScore("frameshift_variant", 0.4);
		makeScore("initiator_codon_variant", 0.3);
		makeScore("intergenic_region", 0.01);
		makeScore("intragenic_variant", 0.001);
		makeScore("intron_variant", 0.05);
		makeScore("missense_variant", 0.9).andMask("protein_altering");
		makeScore("non_coding_transcript_exon_variant", 0.1);
		makeScore("non_coding_transcript_variant", 0.1);
		makeScore("splice_acceptor_variant", 0.5).andMask("protein_altering", "splice");
		makeScore("splice_donor_variant", 0.5).andMask("protein_altering", "splice");
		makeScore("splice_region_variant", 0.5).andMask("protein_altering", "splice");
		makeScore("start_lost", 0.6).andMask("protein_altering");
		makeScore("stop_gained", 0.9).andMask("protein_altering");
		makeScore("stop_lost", 0.6).andMask("protein_altering");
		makeScore("stop_retained_variant", 0.2);
		makeScore("synonymous_variant", 0.1);
		makeScore("upstream_gene_variant", 0.1).andMask("upstream", "updownstream");
	}

	private void dump(SortingCollection<String> sorter, VariantContext ctx) throws Exception {
		if(ctx.getAttributeAsDouble("AF",1.0) >= 0.01) return;

		final List<String> anns = ctx.getAttributeAsStringList(ANN, ".");
		final List<List<String>> predictions = anns.stream().map(PRED -> Arrays.asList(CharSplitter.PIPE.split(PRED)))
				.filter(L -> L.size() > 1).filter(L -> !L.get(1).equals("intergenic_region"))
				.collect(Collectors.toList());

		final String altstr = ctx.getAlternateAllele(0).getDisplayString();
		Double cadd_phred = null;
		Double best_score = null;
		Prediction best_pred = null;
		String best_gene = null;
		// un seul gene par position sinon ca plante

		if (ctx.hasAttribute(CADD_PHRED)) {
			String s = ctx.getAttributeAsString(CADD_PHRED, ".");
			if (s.equals(".") || StringUtils.isBlank(s)) {
				cadd_phred = Double.valueOf(s);
			}
		}

		for (List<String> pred : predictions) {
			if (!pred.get(0).equalsIgnoreCase(altstr))
				continue;
			for (String pred_key : CharSplitter.AMP.split(pred.get(1))) {
				if (pred_key.equals("intergenic_region"))
					continue;
				final Prediction p = scores.getOrDefault(pred_key, null);
				if (p == null)
					throw new IOException(String.join("|", pred));
				if (best_score == null || best_score.compareTo(p.score) < 0) {
					best_score = p.score;
					best_pred = p;
					best_gene = pred.get(3);
				}
			}
		}

		if (best_score != null) {
			best_pred.found = true;

			final StringBuilder sb = new StringBuilder();
			sb.append(ctx.getContig());
			sb.append(" ");
			sb.append(ctx.getStart());
			sb.append(" ");
			sb.append(ctx.getID());
			sb.append(" ");
			sb.append(best_gene);
			sb.append(" ");
			sb.append(best_pred.name);
			sb.append(" ");
			sb.append(best_score);
			sb.append(" ");
			if (cadd_phred == null) {
				sb.append(0.0);
			} else {
				sb.append(cadd_phred.doubleValue());
			}
			sorter.add(sb.toString());
		}

	}

	private static int compareX(String a, String b, int level) {
		String[] t1 = CharSplitter.SPACE.split(a);
		String[] t2 = CharSplitter.SPACE.split(b);
		int i = t1[0].compareTo(t2[0]);// CHROMOSOME
		if (i != 0)
			return i;
		i = t1[3].compareTo(t2[3]);// GENE
		if (i != 0)
			return i;
		if (level == 1)
			return 0;
		return Integer.compare(Integer.parseInt(t1[1]), Integer.parseInt(t2[1]));
	}

	@Override
	public int doWork(List<String> args) {
		try {
			final SortingCollection<String> sorter = SortingCollection.newInstance(String.class, new RowCodec(),
					(A, B) -> compareX(A, B, 2), writingSortingCollection.getMaxRecordsInRam(),
					writingSortingCollection.getTmpPaths());
			sorter.setDestructiveIteration(true);

			try (VCFIterator iter = new VCFIteratorBuilder().open(System.in)) {
				while (iter.hasNext()) {
					final VariantContext vc = iter.next();
					if (vc.getNAlleles() != 2)
						throw new IOException(vc.getContig() + ":" + vc.getStart() + ":" + vc.getAlleles());
					dump(sorter, vc);
				} // end while
			}
			sorter.doneAdding();

			try (CloseableIterator<String> iter0 = sorter.iterator()) {
				try (EqualRangeIterator<String> iter = new EqualRangeIterator<>(iter0, (A, B) -> compareX(A, B, 1))) {

					try (PrintWriter annotOut = IOUtils.openPathForPrintWriter(this.annotationFileOut)) {
						try (PrintWriter setFileOut = IOUtils.openPathForPrintWriter(this.setListFileOut)) {

							while (iter.hasNext()) {
								final List<String> gene_variants = iter.next();

								for (String s : gene_variants) {
									annotOut.println(Arrays.stream(CharSplitter.SPACE.split(s))
											.skip(2L /* chrom pos */)
											.collect(Collectors.joining(" "))
											);
									}
								final String first = gene_variants.get(0);
								final String[] tokens = CharSplitter.SPACE.split(first);
								setFileOut.print(tokens[3]);// gene
								setFileOut.print("\t");
								setFileOut.print(tokens[0]);// contig
								setFileOut.print("\t");
								setFileOut.print((int) gene_variants.stream()
										.mapToInt(it -> Integer.parseInt(CharSplitter.SPACE.split(it)[1])).average()
										.getAsDouble());
								setFileOut.print("\t");
								setFileOut.print(gene_variants.stream().map(it -> CharSplitter.SPACE.split(it)[2])
										.collect(Collectors.joining(",")));
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
					maskOut.print(scores.values().stream().filter(P -> P.masks.contains(mask_name)).map(P -> P.name)
							.collect(Collectors.joining(",")));
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
