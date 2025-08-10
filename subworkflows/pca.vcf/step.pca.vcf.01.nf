/**

Thank you Floriane Simonet for the Help

*/


/**

step_id and step_name are used to get a unique id for multiqc when using this submodule multiple times

*/
if(!params.containsKey("step_id")) throw new IllegalStateException("undefined params.step_id with include+addParams");
if(!params.containsKey("step_name")) throw new IllegalStateException("undefined params.step_name with include+addParams");


include {runOnComplete;dumpParams} from '../../modules/utils/functions.nf'
include {PIHAT01} from '../../subworkflows/pihat/pihat.01.nf' 


def whatisapca = "<cite><b>PCA</b> is a statistical technique for reducing the dimensionality of a dataset. This is accomplished by linearly transforming the data into a new coordinate system where (most of) the variation in the data can be described with fewer dimensions than the initial data.</cite>"


workflow ACP_VCF_STEP {
	take:
		genome_ch
		vcf
		sample2collection
		blacklisted_bed
		exclude_samples
	main:
		pihat_ch = PIHAT01(genome_ch, vcf, exclude_samples , blacklisted_bed)
		cluster_ch = PLINK_CLUSTER( pihat_ch.genome_bcf, pihat_ch.plink_genome)
		
		assoc_ch = PLINK_ASSOC(genome_ch , pihat_ch.genome_bcf, cluster_ch.output)

		plot_assoc_ch = PLOT_ASSOC(genome_ch, assoc_ch.assoc.flatten())

		CLEANUP_VCF(vcf, assoc_ch.variants_to_remove)
	

		headers = cluster_ch.header.splitCsv(header:true,sep:'\t')
		
		plot_ch = PLOT_IT(
			headers.combine(headers).
			filter{T->T[0].label.compareTo(T[1].label)<0}.
			map{T->[
				label1:T[0].label,
				label2:T[1].label,
				]}.
			combine(cluster_ch.output).
			map{T->T[0].plus(clusters:T[1])},
			sample2collection
			)


		to_multiqc = plot_ch.output.mix(plot_assoc_ch.output.flatten()).mix(pihat_ch.to_multiqc)

	emit:
		multiqc = to_multiqc

	}

process PLINK_CLUSTER {
	tag "${genome_bcf.name} ${genome_plink.name}"
	conda "${moduleDir}/../../conda/bioinfo.01.yml"

	input:
		path(genome_bcf)
		path(genome_plink)
	output:
		path("cluster.mds"),emit:output
		path("header.tsv"),emit:header
	script:
		def num_components = 3
	"""
	hostname 2>&1
	mkdir -p TMP

	plink --bcf '${genome_bcf}' \\
		--double-id \\
		--read-genome '${genome_plink}' \\
		--mds-plot ${num_components} \\
		--cluster \\
		--out TMP/cluster

	# reformat mds, normalize spaces, replace with tab
	tr -s " " < TMP/cluster.mds | sed -e 's/^[ ]*//' -e 's/[ ]*\$//' | tr " " "\t" > cluster.mds

	# column 1 is FID
	# column 2 is IID
	# remaining are 

	head -n1 cluster.mds |\\
		tr "\t" "\\n" |\\
		awk -F '\t' 'BEGIN {printf("label\\n");} (\$1 ~ /^C[0-9]*\$/ ) {printf("%s\\n",\$1)}' > header.tsv
	"""
	}

process PLINK_ASSOC {
	tag "${genome_bcf.name} ${cluster_mds.name}"
	conda "${moduleDir}/../../conda/bioinfo.01.yml"
	input:
		tuple path(fasta),path(fai),path(dict)
		path(genome_bcf)
		path(cluster_mds)
	output:
		path("PHENO.*.qassoc"),emit:assoc
		path("variants_to_remove.bcf"),emit:variants_to_remove
		path("variants_to_remove.bcf.csi"),emit:variants_to_remove_csi
	script:
		def treshold = 5E-8
	"""
	hostname 2>&1
	mkdir -p TMP

	# create pheno file https://www.cog-genomics.org/plink/1.9/input#pheno
	awk '(NR==1) {split(\$0,header);} {X=0;for(i=1;i<=NF;i++) {if(header[i] !="SOL") {printf("%s%s",(X==0?"":"\t"),\$i);X=1;}} printf("\\n");}' '${cluster_mds}' > TMP/pheno.tsv

	plink --bcf  '${genome_bcf}' \\
		--double-id \\
		--allow-no-sex \\
		--pheno TMP/pheno.tsv \\
		--all-pheno \\
		--assoc \\
		--out PHENO
	
	# list position CHROM POS REF ALT to remove
	echo "##fileformat=VCFv4.2" > TMP/jeter.vcf
	echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> TMP/jeter.vcf

	LC_ALL=C awk '(\$NF!="P" && \$NF < ${treshold}) {print \$2}'  PHENO.*.qassoc  |\\
		awk -F ':' '{printf("%s\t%s\t.\t%s\t%s\t.\t.\t.\\n",\$1,\$2,\$3,\$4);}' |\
		sed 's/^chr//' >> TMP/jeter.vcf
	
	jvarkit vcfsetdict -R '${fasta}'  --onNotFound SKIP TMP/jeter.vcf |\
		bcftools sort -T TMP/tmp -O b -o variants_to_remove.bcf
	bcftools index -f variants_to_remove.bcf
	"""
	}

process CLEANUP_VCF {
	input:
		path(vcf)
		path(variants_to_remove)
	output:
	script:
	"""
	"""
}

process PLOT_ASSOC {
tag "${assoc.name}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	tuple path(fasta),path(fai),path(dict)
	path(assoc)
output:
	path("${prefix}.multiqc.*"),emit:output
script:
	prefix = "plotassoc"
	def title = assoc.name.replace('.','_')
	def nhead = 10
"""
hostname 1>&2
mkdir -p TMP

JD1=`which jvarkit`
echo "JD1=\${JD1}" 1>&2
# directory of jvarkit
JD2=`dirname "\${JD1}"`
# find the jar itself
JVARKIT_JAR=`find "\${JD2}/../.." -type f -name "jvarkit.jar" -print -quit`
JVARKIT_DIST=`dirname "\${JVARKIT_JAR}"`


cat << __EOF__ > TMP/Minikit.java
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
private final static double MIN_P = ${params.plot_min_pvalue?:"1E-100"};
private int doWork() {
    final String REF="${fasta}";
    final String input = "${assoc}";
    try {
        Path path = Paths.get(input);
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
__EOF__

javac -d TMP -cp \${JVARKIT_DIST}/jvarkit.jar TMP/Minikit.java
java -Djava.awt.headless=true -cp \${JVARKIT_DIST}/jvarkit.jar:TMP Minikit


cat << __EOF__ > TMP/jeter.html
<!--
parent_id: pihat_section
parent_name: "PCA ${params.step_name}"
parent_description: "${whatisapca}"
id: '${title}_table'
section_name: '${title} table'
description: '${assoc} :  ${nhead} first lines. Number of variants per Chromosome.'
-->
__EOF__


tr -s " " < "${assoc}" | LC_ALL=C sort -T TMP -t ' ' -k10,10g  | head -n ${nhead} |\\
	awk 'BEGIN{printf("<table class=\\"table\\">");} {t=(NR==1?"th":"td");print("<tr>");for(i=1;i<=NF;i++) printf("<%s>%s</%s>",t,\$i,t);printf("</tr>\\n");} END{printf("</table>\\n");}'  >> TMP/jeter.html

#
# Floriane suggested to also add the number of variants per chromosome
#
awk '(\$2!="SNP") {print \$2}' '${assoc}' |\\
	cut -d ':' -f 1 | sort -T TMP | uniq -c |\\
	awk '{printf("%s\t%s\\n",\$2,\$1);}' |\\
	sort -T TMP -t '\t' -k1,1 |\\
	join -t '\t' -1 1 -2 1 - <(sort -t '\t' -k1,1 '${fai}' ) |\\
	awk 'BEGIN{printf("<table class=\\"table\\"><thead><caption>Number of variants per contig</caption><tr><th>CHROM</th><th>Length</th><th>VARIANTS</th><th>Variants per base</th></tr></thead><tbody>\\n");} {printf("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>\\n",\$1,\$3,\$2,\$2/\$3);} END {printf("</tbody></table>\\n");}' >> TMP/jeter.html



mv TMP/jeter.html "${prefix}.multiqc.${title}.table_mqc.html"


##
## create MULTIQC CONFIG
##
cat << EOF > "${prefix}.multiqc.${title}.multiqc_config.yaml"
custom_data:
  pca_${title}:
    parent_id: pihat_section
    parent_name: "PCA ${params.step_name}"
    parent_description: "${whatisapca}"
    section_name: "PCA: ${assoc.name}"
    description: "PCA: ${assoc.name}"
sp:
  pca_${title}:
    fn: "${prefix}.multiqc.${title}.png"
ignore_images: false
EOF

"""
}

process PLOT_IT {
	tag "${row.label1} vs ${row.label2}"
	conda "${moduleDir}/../../conda/bioinfo.01.yml"

	input:
		val(row)
                path(sample2collection)
	output:
		path("${prefix}.multiqc.*"),emit:output
	script:
		prefix = "plot"
		def title = row.label1+"_"+row.label2
	"""
	hostname 2>&1
	mkdir -p TMP

echo "IID\tcollection" > TMP/sn2col.tsv

if ${!sample2collection.name.equals("NO_FILE")} ; then

	#  plink asshole for naming samples...
	awk -F '\t' '/^#/ {next} {printf("%s_%s_%s_%s\t%s\\n",\$1,\$1,\$1,\$1,\$2);}' '${sample2collection}' |\\
		sort | uniq >> TMP/sn2col.tsv

fi


cat << 'EOF' > TMP/jeter.R
# load sample/collection
sn2col <- read.table("TMP/sn2col.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

data <- read.table("${row.clusters}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

mergedata <- merge(data, sn2col, by = "IID", all.x = TRUE)

mergedata\$collection[is.na(mergedata\$collection)] <- "OTHER"

head(mergedata)

COLORS <- rainbow(length(unique(mergedata\$collection)),alpha = 0.75)

# Sélectionner les colonnes C*
colX <- mergedata\$${row.label1}
colY <- mergedata\$${row.label2}
pigments <- COLORS[as.factor(mergedata\$collection)]

# Créer le nuage de points
png("${prefix}.multiqc.${title}.pca.png") 
plot(colX, colY, main = "${row.label1} vs ${row.label2}",
	sub= "${row.clusters}",
	xlab = "${row.label1}",
	ylab = "${row.label2}",
	pch = 16,
	col = pigments
	)

if (${!sample2collection.name.equals("NO_FILE")?"TRUE":"FALSE"})  {

	legend("topright", legend = levels(as.factor(mergedata\$collection)), col = COLORS, pch = 16, title = "Collection", cex = 0.7)
}


dev.off()  # Fermeture du fichier de sortie
EOF

R --vanilla < TMP/jeter.R

##
## create MULTIQC CONFIG
##
cat << EOF > "${prefix}.multiqc.config.yaml"
custom_data:
  pca_${title}:
    parent_id: pihat_section
    parent_name: "PCA ${params.step_name}"
    parent_description: "${whatisapca}"
    section_name: "PCA: ${row.label1} ${row.label2}"
    description: "PCA: ${row.label1} vs ${row.label2}"
sp:
  pca_${title}:
    fn: "${prefix}.multiqc.${title}.pca.png"
ignore_images: false
EOF

	"""
}
