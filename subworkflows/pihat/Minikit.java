import java.awt.Color;
import java.awt.Font;
import java.awt.AlphaComposite;
import java.awt.Composite;
import java.awt.Graphics2D;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.imageio.ImageIO;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;


public class Minikit {
private final static double MIN_P =__MIN_P_VALUE__;
private int doWork() {
    final String REF="${fasta}";
    final String input = "${assoc}";
    try {
        Path path = Paths.get("__INPUT__");
        final SAMSequenceDictionary dict0 = SAMSequenceDictionaryExtractor.extractDictionary(Paths.get(REF));
        final Pattern regex = Pattern.compile("(chr)?[XY0-9]+");
        final SAMSequenceDictionary dict = new SAMSequenceDictionary(dict0.getSequences().stream().
                filter(SR->regex.matcher(SR.getSequenceName()).matches()).
                collect(Collectors.toList()));
        final Pattern ws = Pattern.compile("\\\\s+");
        final Pattern colon = Pattern.compile(":");
        final double p_treshold = -Math.log10(5E-8);
        double max_p = p_treshold;
        double min_p = 0.0;
        try(BufferedReader br = Files.newBufferedReader(path)) {
            for(;;) {
                final String line = br.readLine();
                if(line==null) break;
                final String[] tokens = ws.split(line.trim());
                if(tokens[tokens.length-1].equals("P")) continue;
                double p = Double.parseDouble(tokens[tokens.length-1]);
		if(p< MIN_P) p = MIN_P;
                p = - Math.log10(p);
                max_p = Math.max(max_p, p);
                min_p = Math.min(min_p, p);
                }
            max_p = max_p + (max_p - min_p)*0.1;
            }
        int width = 1000;
        int height = 300;
        int left_margin = 0;
        int bottom_margin = 0;
        final double genome_len = dict.getReferenceLength();
        final BufferedImage img = new BufferedImage(width+left_margin+1, height+bottom_margin+1, BufferedImage.TYPE_INT_RGB);
        Graphics2D g=(Graphics2D)img.getGraphics();
        g.setColor(Color.WHITE);
        g.fillRect(0,0,img.getWidth(),img.getHeight());



        g.setColor(Color.DARK_GRAY);
        long pos=0;
        for(int tid=0;tid < dict.size();++tid) {
            final SAMSequenceRecord ssr = dict.getSequence(tid);
            final double x = left_margin + (pos/genome_len)*width;
            g.setColor(Color.DARK_GRAY);
	    g.setFont(new Font(g.getFont().getName(),Font.PLAIN,10));
	    g.drawString(ssr.getSequenceName(),(int)x+3,10);
            g.setColor(Color.DARK_GRAY);
            g.draw(new Line2D.Double(x,0,x,height));
            pos+=ssr.getSequenceLength();
            }
        
        g.setColor(Color.RED);
        final double y_treshold = height - ((p_treshold-min_p)/(max_p-min_p))*height;
        g.draw(new Line2D.Double(left_margin,y_treshold,left_margin+width,y_treshold));
        final double radius=1.5;
        
        try(BufferedReader br = Files.newBufferedReader(path)) {
            for(;;) {
                String line = br.readLine();
                if(line==null) break;
                String[] tokens = ws.split(line.trim());
                if(tokens[tokens.length-1].equals("P")) continue;
                double p = Double.parseDouble(tokens[tokens.length-1]);
		boolean low_pvalue = false;
		if( p < MIN_P) {
			low_pvalue = true;
			p = MIN_P;
			}
                p = - Math.log10(p);
                tokens = colon.split(tokens[1]);
                final String contig = tokens[0];
                pos = 0;
                int tid=0;
                Color col = Color.BLACK;
                for(tid=0;tid < dict.size();++tid) {
                    final SAMSequenceRecord ssr = dict.getSequence(tid);
                    if(ssr.getSequenceName().equals(contig)) {
                        if(contig.equals("X") || contig.equals("chrX")) {
                            col = Color.BLUE;
                            }
                        else if(contig.equals("Y") || contig.equals("chrY")) {
                            col = Color.PINK;
                            }
                        else
                            {
                            col = (tid%2==0?Color.GREEN:Color.MAGENTA);
                            }
                        pos+= Integer.parseInt(tokens[1]);
                        break;
                        }
                    pos+=ssr.getSequenceLength();
                    }
                if(tid==dict.size()) {
			continue;
			}
                double x = left_margin + (pos/genome_len)*width;
                double y = height - ((p-min_p)/(max_p-min_p))*height;
                g.setColor(col);
		final Composite oldcomposite = g.getComposite();
		if(low_pvalue) {
			final int radius2=5;
			g.setColor(Color.MAGENTA);
			g.fill(new Ellipse2D.Double(x-radius2,0-radius2,radius2*2,radius2*2));
			}
                else if(p> p_treshold) {
                    g.fill(new Ellipse2D.Double(x-radius,y-radius,radius*2,radius*2));
                    }
                else {
		    g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER,0.5f));
                    g.draw(new Line2D.Double(x-radius,y,x+radius,y));
                    g.draw(new Line2D.Double(x,y-radius,x,y+radius));
                    }
		g.setComposite(oldcomposite);
                }
            }
        
	g.setColor(Color.BLACK);        
        g.draw(new Rectangle2D.Double(left_margin,0,width,height));
        g.dispose();
        ImageIO.write(img, "PNG", new File("${prefix}.multiqc.${title}.png"));
        return 0;
        }
    catch(final Throwable err) {
        err.printStackTrace();
        return -1;
        }
    }
    
public static void main(final String[] args)
    {
    int ret= new Minikit().doWork();
    System.exit(ret);
    }
}