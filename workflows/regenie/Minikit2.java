import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class Minikit extends Launcher {
	@Parameter(names={"-o","--output"},description= "output dir",required = true)
	private Path outputDir = null;
	@Parameter(names={"-R"},description= INDEXED_FASTA_REFERENCE_DESCRIPTION,required = true)
	private Path reference = null;

	private static class Condition {
		Path path;
		PrintWriter w;
		String annot="";
		String freq_name;
		String test_name;
		int count=0;
		}
	@Override
	public int doWork(List<String> args) {
		try {
			final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(reference);
			if(dict==null) throw new IllegalArgumentException();
			final String kX  = dict.getSequences().stream().map(SSR->SSR.getSequenceName()).filter(S->S.matches("(chr)?X")).findFirst().orElse(null);
			final String kY  = dict.getSequences().stream().map(SSR->SSR.getSequenceName()).filter(S->S.matches("(chr)?Y")).findFirst().orElse(null);
			ContigNameConverter convert = ContigNameConverter.fromOneDictionary(dict);
			
			String input = oneFileOrNull(args);
			final List<Condition> conditions = new ArrayList<>();
			try(BufferedReader br = super.openBufferedReader(input)) {
				String line = br.readLine();
				if(line==null) throw new IOException("cannot read first line");
				final String header_line = line;
				final FileHeader header = new FileHeader(line, CharSplitter.SPACE);
				final int col_test = header.getColumnIndex("TEST");
				final int col_id = header.getColumnIndex("ID");
				final int col_chrom = header.getColumnIndex("CHROM");
				final int col_extra = header.getColumnIndex("EXTRA");
				final int col_log10P = header.getColumnIndex("LOG10P");
				while((line=br.readLine())!=null) {
					final List<String> tokens = header.split(line);
					if(tokens.get(col_log10P).equals("NA")) continue;
					if(tokens.get(col_extra).equals("TEST_FAIL")) continue;
					final String test_name = tokens.get(col_test);
					final String id =  tokens.get(col_id);
					final String freq_name;
					String s3;
					if(id.endsWith(".singleton")) {
						freq_name = "singleton";
						s3 = id.substring(0,id.length()-10);
						}
					else
						{
						int x = id.lastIndexOf("0.");
						if(x==-1) throw new IllegalArgumentException("no 0. in "+id);
						if(id.charAt(x-1)!='.') throw new IllegalArgumentException("no . in "+id);
						s3 =  id.substring(0,x-1);//-1 for dot before the number
						freq_name = id.substring(x);
						//check double
						Double.parseDouble(freq_name);
						}
					int dot2 = s3.lastIndexOf('.');
					if(dot2==-1) throw new IllegalArgumentException("no '.' in "+s3);
					final String annot= s3.substring(dot2+1);
					if(annot.isEmpty())  throw new IllegalArgumentException("empty for /"+s3+"/"+annot+"/"+id);

					Condition cond = conditions.stream().filter(C->C.annot.equals(annot) && C.test_name.equals(test_name) && C.freq_name.equals(freq_name)).findFirst().orElse(null);
					if(cond==null) {
						cond = new Condition();
						cond.path =outputDir.resolve(String.join("_",freq_name,test_name,annot)+".regenie");
						cond.test_name = test_name;
						cond.freq_name = freq_name;
						cond.annot = annot;
						cond.w= IOUtils.openPathForPrintWriter(cond.path);
						cond.w.println(header_line);
						conditions.add(cond);
						}
					if(kX!=null && tokens.get(col_chrom).equals("23")) {
						tokens.set(col_chrom,kX);
						}
					else if(kY!=null && tokens.get(col_chrom).equals("24")) {
						tokens.set(col_chrom,kY);
						}
					else
						{
						final String c = convert.apply(tokens.get(col_chrom));
						if(StringUtils.isBlank(c)) {
							System.err.println("unknown contig "+tokens.get(col_chrom));
							continue;
							}
						tokens.set(col_chrom,c);
						}
					cond.w.println(String.join(" ",tokens));
					cond.count++;
					}
				
				try(PrintWriter manifest = super.openPathOrStdoutAsPrintWriter(this.outputDir.resolve("manifest.tsv"))) {
					manifest.println("filename\tfreq\ttest_name\tannot");
					for(Condition cond : conditions) {
						if(cond.count>0) {
							manifest.print(cond.path);
							manifest.print("\t");
							manifest.print(cond.freq_name);
							manifest.print("\t");
							manifest.print(cond.test_name);
							manifest.print("\t");
							manifest.print(cond.annot);
							manifest.println();
							}
						cond.w.flush();
						cond.w.close();
						}
					manifest.flush();
					}
				}
			return 0;
			}
		catch(Throwable err) {
			err.printStackTrace();
			return -1;
			}
		}
	
	public static void main(final String[] args) {
		new Minikit().instanceMainWithExit(args);
	}
}
