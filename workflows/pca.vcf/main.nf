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


		to_multiqc = plot_ch.multiqc.mix(plot_ch.multiqc_yaml)


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
	script:
		def num_components = 3
	"""
	hostname 2>&1
	${moduleLoad("plink bcftools r/3.6.3")}
	mkdir -p TMP

	plink --bcf '${genome_bcf}' \\
		--double-id \\
		--read-genome '${genome_plink}' \\
		--mds-plot ${num_components} \\
		--cluster \\
		--out TMP/cluster


	awk '(NR==1) {split(\$0,header);next;} {X=0; for(i=1;i<=NF;i++) {if(header[i] ~ /^C[0-9]+\$/) {printf("%s%s",(X==0?"":"\t"),\$i);X=1;}} printf("\\n");}' TMP/cluster.mds > cluster.tsv

	head -n1 cluster.tsv |\\
		tr "\t" "\\n" |\\
		awk -F '\t' 'BEGIN {printf("column\tlabel\\n");} {printf("%d\tC%d\\n",NR,NR)}' > header.tsv
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
