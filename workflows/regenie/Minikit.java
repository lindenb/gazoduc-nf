
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class Minikit {
	private Map<String,Double> scores= new HashMap<>();
	private final Pattern amp = Pattern.compile("[&]"); 
	public Minikit() {
		scores.put("3_prime_UTR_variant", 0.1);
		scores.put("5_prime_UTR_premature_start_codon_gain_variant", 0.2);
		scores.put("5_prime_UTR_truncation", 0.5);
		scores.put("5_prime_UTR_variant", 0.2);
		scores.put("bidirectional_gene_fusion", 1.0);
		scores.put("conservative_inframe_deletion", 0.1);
		scores.put("conservative_inframe_insertion", 0.3);
		scores.put("disruptive_inframe_deletion", 0.2);
		scores.put("disruptive_inframe_insertion", 0.2);
		scores.put("downstream_gene_variant", 0.1);
		scores.put("exon_loss_variant", 1.0);
		scores.put("frameshift_variant", 0.4);
		scores.put("initiator_codon_variant", 0.3);
		scores.put("intergenic_region", 0.01);
		scores.put("intragenic_variant", 0.001);
		scores.put("intron_variant", 0.05);
		scores.put("missense_variant", 0.9);
		scores.put("non_coding_transcript_exon_variant", 0.1);
		scores.put("non_coding_transcript_variant", 0.1);
		scores.put("splice_acceptor_variant", 0.5);
		scores.put("splice_donor_variant", 0.5);
		scores.put("splice_region_variant", 0.5);
		scores.put("start_lost", 0.6);
		scores.put("stop_gained", 0.9);
		scores.put("stop_lost", 0.6);
		scores.put("stop_retained_variant", 0.2);
		scores.put("synonymous_variant", 0.1);
		scores.put("upstream_gene_variant", 0.1);
		}
	
	private void dump(String contig,String pos, String id, String ref, String alt,List<List<String>> predictions)  throws Exception  {
		
			Double best_score=null;
			String best_pred=null;
			String best_gene = null;
			// un seul gene par position sinon ca plante
			
			for(List<String> pred : predictions) {
				if(!pred.get(0).equalsIgnoreCase(alt)) continue;
				for(String pred_key: this.amp.split(pred.get(1))) {
					Double v = scores.getOrDefault(pred_key, null);
					if(v==null) throw new IOException(String.join("|",pred));
					if(best_score==null || best_score.compareTo(v)<0) {
						best_score=v;
						best_pred = pred_key;
						best_gene = pred.get(3);
						}
					}
	 			}
			
			if(best_score!=null) {
				System.out.print(contig);
                                System.out.print(" ");
				System.out.print(pos);
                                System.out.print(" ");
				System.out.print(id);
				System.out.print(" ");
				System.out.print(best_gene);
				System.out.print(" ");
				System.out.print(best_pred);
				System.out.print(" ");
				System.out.print(best_score);
				System.out.println();
				}
			
		}
	private void run()  throws Exception {
		final Pattern semicolon = Pattern.compile("[;]"); 
		final Pattern tab = Pattern.compile("[\t]"); 
		final Pattern comma = Pattern.compile("[,]");
		final Pattern pipe = Pattern.compile("[\\|]");
		
		try(BufferedReader br= new BufferedReader(new InputStreamReader(System.in))) {
			String line;
			while((line=br.readLine())!=null) {
				if(line.startsWith("#")) continue;
				String[] tokens = tab.split(line);
				if(tokens[4].contains(","))  throw new IOException(tokens[4]);
				String ann = Arrays.stream(semicolon.split(tokens[7])).filter(S->S.startsWith("ANN=")).findFirst().orElse(null);
				if(ann==null) continue;
				ann=ann.substring(4);
				final List<List<String>> predictions= Arrays.stream(comma.split(ann)).
						map(PRED->Arrays.asList(pipe.split(PRED))).
						filter(L->!L.get(1).equals("intergenic_region")).
						collect(Collectors.toList());
				dump(tokens[0],tokens[1],tokens[2],tokens[3],tokens[4],predictions);
				}
			}
		}

	public static void main(final String[] args)
		{
		try {
				new Minikit().run();
			}
		catch(Throwable err) {
			err.printStackTrace();
			System.exit(-1);
			}
		}
	}
