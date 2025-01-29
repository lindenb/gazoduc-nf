import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;

import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Locatable;


public class Minikit extends Launcher {
	@Parameter(names={"-o","--out"},description=ArchiveFactory.OPT_DESC,required=true)
	private Path outputFile= null;
	@Parameter(names={"-w"},description="window size"+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.LongStringConverter.class,splitter=NoSplitter.class)
	private long window_size=-1L;
	@Parameter(names={"-s"},description="window shift"+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.LongStringConverter.class,splitter=NoSplitter.class)
	private long window_shift=-1L;

	
	private int id_generator =0;
	

	private void dump(
			final ArchiveFactory archiveFactory,
			final List<Locatable> src
			) throws IOException {
		if(src.isEmpty()) return;
		final int chromStart = src.stream().mapToInt(R->R.getStart()-1).min().orElse(0);
		final int chromEnd = src.stream().mapToInt(R->R.getEnd()).max().orElse(0);
		final String prefix= src.get(0).getContig()+"_"+(chromStart+1)+"_"+chromEnd;
		final String filename = String.format("%s.%09d.bed",prefix, ++id_generator);

		
		try(final PrintWriter pw =  new PrintWriter(archiveFactory.openWriter(filename))) {
			for(final Locatable r: src) {
				pw.print(r.getContig());
				pw.print("\t");
				pw.print(r.getStart()-1);//convert to bed if needed
				pw.print("\t");
				pw.print(r.getEnd());
				pw.println();
				}
			pw.flush();
			}


		}
	private void apply_cluster(
			final ArchiveFactory archiveFactory,
			final List<Locatable> src
			) throws IOException
		{
		Collections.sort(src,(A,B)->{
			int i = Integer.compare(A.getStart(), B.getStart());
			if(i!=0) return i;
			return  Integer.compare(A.getEnd(), B.getEnd());
			});
		
		int prev_start=-1;
		for(int x=0;x< src.size();x++) {
			final Locatable loca = src.get(x);
			if(prev_start>-1 && prev_start+this.window_shift > loca.getStart()) continue;
			final List<Locatable> cluster=new ArrayList<>();
			cluster.add(loca);
			prev_start = loca.getStart();
			for(int y=x+1;y< src.size();y++) {
				final Locatable locb = src.get(y);
				if(CoordMath.getLength(loca.getStart(), locb.getEnd())> this.window_size) break;
				cluster.add(locb);
				}
			dump(archiveFactory,cluster);
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		if(this.window_size<0) return -1;
		if(this.window_shift<0) return -1;
		try
			{
			try(ArchiveFactory archiveFactory = ArchiveFactory.open(this.outputFile)) {
				final String input = oneFileOrNull(args);
				final Map<String,List<Locatable>> contig2lines = new HashMap<>();
				try(BedLineReader br = new BedLineReader(super.openBufferedReader(input),input)) {
					br.stream().
						map(R->new SimpleInterval(R)).
						forEach(RGN->{
							List<Locatable> L = contig2lines.get(RGN.getContig());
							if(L==null) {
								L = new ArrayList<>();
								contig2lines.put(RGN.getContig(),L);
								}
							L.add(RGN);
						});
					}
				for(List<Locatable> losc: contig2lines.values()) {
					apply_cluster(archiveFactory, losc);
					}
				}
		
			return 0;
			}
		catch(final Throwable err)
			{
			err.printStackTrace();
			return -1;
			}
		}
	

	public static void main(final String[] args)
		{
		new Minikit().instanceMainWithExit(args);
		}
	}
