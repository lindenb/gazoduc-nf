/**

Thank you Floriane Simonet for the Help

*/

include {moduleLoad;runOnComplete;dumpParams} from '../../modules/utils/functions.nf'
include {PIHAT01} from '../../subworkflows/pihat/pihat.01.nf'
include {MULTIQC_01} from '../../modules/multiqc/multiqc.01.nf'
include {PARAMS_MULTIQC} from '../../modules/utils/params.multiqc.nf'


if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}


runOnComplete(workflow)

workflow {
	ACP_VCF([:], params.genomeId, file(params.vcf))
	}

workflow ACP_VCF {
	take:
		meta
		genomeId
		vcf
	main:
		pihat_ch = PIHAT01(genomeId, vcf, Channel.fromPath("NO_FILE"))
		cluster_ch = PLINK_CLUSTER( pihat_ch.genome_bcf, pihat_ch.plink_genome)
		
		assoc_ch = PLINK_ASSOC( genomeId, pihat_ch.genome_bcf, cluster_ch.mds)

		plot_assoc_ch = PLOT_ASSOC(genomeId, assoc_ch.assoc.flatten())

		CLEANUP_VCF(vcf, assoc_ch.variants_to_remove)
	

		headers = cluster_ch.header.splitCsv(header:true,sep:'\t')
		
		plot_ch = PLOT_IT(headers.combine(headers).
			filter{T->T[0].label.compareTo(T[1].label)<0}.
			map{T->[
				column1:T[0].column,
				label1:T[0].label,
				column2:T[1].column,
				label2:T[1].label,
				]}.
			combine(cluster_ch.output).
			map{T->T[0].plus(clusters:T[1])})


		to_multiqc = plot_ch.multiqc.mix(plot_ch.multiqc_yaml).mix(plot_assoc_ch.output)


                params4multiqc_ch = PARAMS_MULTIQC([:])
                mqc_ch = MULTIQC_01(
                                ["title":"${params.prefix}PCA"],
                                to_multiqc.
				concat(params4multiqc_ch.output).
				concat(pihat_ch.to_multiqc).
				collect()
                                );

	}

process PLINK_CLUSTER {
	tag "${genome_bcf.name} ${genome_plink.name}"
	input:
		path(genome_bcf)
		path(genome_plink)
	output:
		path("cluster.tsv"),emit:output
		path("header.tsv"),emit:header
		path("cluster.mds"),emit:mds
	script:
		def num_components = 3
	"""
	hostname 2>&1
	${moduleLoad("plink")}
	mkdir -p TMP

	plink --bcf '${genome_bcf}' \\
		--double-id \\
		--read-genome '${genome_plink}' \\
		--mds-plot ${num_components} \\
		--cluster \\
		--out TMP/cluster

	mv TMP/cluster.mds ./

	awk '(NR==1) {split(\$0,header);next;} {X=0; for(i=1;i<=NF;i++) {if(header[i] ~ /^C[0-9]+\$/) {printf("%s%s",(X==0?"":"\t"),\$i);X=1;}} printf("\\n");}' cluster.mds > cluster.tsv

	head -n1 cluster.tsv |\\
		tr "\t" "\\n" |\\
		awk -F '\t' 'BEGIN {printf("column\tlabel\\n");} {printf("%d\tC%d\\n",NR,NR)}' > header.tsv
	"""
	}

process PLINK_ASSOC {
	tag "${genome_bcf.name} ${cluster_mds.name}"
	input:
		val(genomeId)
		path(genome_bcf)
		path(cluster_mds)
	output:
		path("PHENO.*.qassoc"),emit:assoc
		path("variants_to_remove.bcf"),emit:variants_to_remove
		path("variants_to_remove.bcf.csi"),emit:variants_to_remove_csi
	script:
		def treshold = 5E-8
		def reference = params.genomes[genomeId].fasta
	"""
	hostname 2>&1
	${moduleLoad("plink bcftools jvarkit")}
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
	
	java -jar \${JVARKIT_DIST}/jvarkit.jar vcfsetdict -R '${reference}'  --onNotFound SKIP TMP/jeter.vcf |\
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
input:
	val(genomeId)
	path(assoc)
output:
	path("multiqc_config.yaml"),emit:output
script:
	def reference = params.genomes[genomeId].fasta
	def title = assoc.name.replace('.','_')
"""
hostname 1>&2
${moduleLoad("jvarkit")}
mkdir -p TMP

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

private int doWork() {
    final String REF="${reference}";
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
                if(p> p_treshold) {
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
        ImageIO.write(img, "PNG", new File("${title}.png"));
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
java -cp \${JVARKIT_DIST}/jvarkit.jar:TMP Minikit


##
## create MULTIQC CONFIG
##
cat << EOF > multiqc_config.yaml
custom_data:
  pca_${title}:
    parent_id: pihat_section
    parent_name: "PCA"
    parent_description: "PCA"
    section_name: "PCA: ${assoc.name}"
    description: "PCA: ${assoc.name}"
sp:
  pca_${title}:
    fn: "\${PWD}/${title}.png"
ignore_images: false
EOF

"""
}

process PLOT_IT {
	tag "${row.label1} vs ${row.label2}"
	input:
		val(row)
	output:
		path("*.png"),emit:multiqc
	        path("multiqc_config.yaml"),emit:multiqc_yaml
	script:

		def title = row.label1+"_"+row.label2
	"""
	hostname 2>&1
	${moduleLoad("r/3.6.3")}
	mkdir -p TMP

cat << 'EOF' > TMP/jeter.R
data <- read.table("${row.clusters}", header = FALSE, sep = "\t")

# Sélectionner les colonnes 3 et 5
colX <- data[, ${row.column1}]
colY <- data[, ${row.column2}]

# Créer le nuage de points
png("${title}.pca.png") 
plot(colX, colY, main = "${row.label1} vs ${row.label2}",
	sub= "${row.clusters}",
	xlab = "${row.label1}",
	ylab = "${row.label2}",
	pch = 16,
	col = "blue"
	)
dev.off()  # Fermeture du fichier de sortie
EOF

R --vanilla < TMP/jeter.R

##
## create MULTIQC CONFIG
##
cat << EOF > multiqc_config.yaml
custom_data:
  pca_${title}:
    parent_id: pihat_section
    parent_name: "PCA"
    parent_description: "PCA"
    section_name: "PCA: ${row.label1} ${row.label2}"
    description: "PCA: ${row.label1} ${row.label2}"
sp:
  pca_${title}:
    fn: "${title}.pca.png"
ignore_images: false
EOF

	"""
}
