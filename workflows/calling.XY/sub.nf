workflow MAKE_STATS {
	take:
		fasta
		fai
		gtf
		samplesheet
		bcf_files
		bed_and_sex_ctg_ch
	main:
		sex_ch = bed_and_sex_ctg_ch.map{it[1]}.unique()

		ch1 = EXTRACT_SAMPLES( samplesheet, sex_ch )

		

	
		ch2 = BCFTOOLS_STATS(fasta, fai, gtf, bcf_files,
			bed_and_sex_ctg_ch.
				combine(ch1.output).
				filter{it[1].equals(it[3])}.
				map{[it[1],it[0],it[2],it[4],it[5]]}
			)
		
		ch3 = MULTIQC_MALE_VS_FEMALE(ch2.directory.filter{it[0].contains("M") && it[0].contains("F")});
		ch4 = ch2.output.mix(ch3.output)
		
	emit:
		output = ch4		
	}


process EXTRACT_SAMPLES {
label "process_single"
tag "${sex}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	path(samplesheet)
	val(sex)
output:
	tuple val(sex),path("samples.${sex}.txt"),path("sample2sex.${sex}.tsv"),emit:output
script:
	def sexm = sex.contains("M")?"male":"_ignore"
	def sexf = sex.contains("F")?"female":"_ignore"
"""
mkdir -p TMP

cat << '__EOF__' > TMP/jeter.R
data <- read.table("${samplesheet}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
filtered_data <- data[data\$sex == "${sexm}" | data\$sex == "${sexf}", ]
filtered_data <- filtered_data[, c("sample", "sex")]
write.table(filtered_data, file = "sample2sex.${sex}.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
__EOF__

R --vanilla --no-save < TMP/jeter.R
cut -f1 sample2sex.${sex}.tsv | sort | uniq > samples.${sex}.txt
"""
}

process BCFTOOLS_STATS {
label "process_single"
tag "sex:${sex} bed:${bed.name} samples:${samples.name} ctg:${contig}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	path(fasta)
	path(fai)
	path(gtf)
	path(bcf_files)
	tuple val(sex),path(bed),val(contig),path(samples),path(sample2sex)
output:
	path("multiqc.${sex}.${bed.name}.${contig}.zip"),emit:output
	tuple val(sex),val("${bed.name}"),val(contig),path("multiqc.${sex}.${bed.name}.${contig}"),path(sample2sex),emit:directory
script:
	def vcf= bcf_files.find{it.name.endsWith("f")}
"""
mkdir -p TMP
set -o pipefail
export LC_ALL=en_US.utf8

gunzip -c "${gtf}" |\\
	awk -F '\t' '(\$3=="exon" && \$1 ~ /^(chr)?[${contig}]\$/)' |\\
	cut -f 1,4,5 |\\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\\
	bgzip > TMP/exons.txt.gz




if ${bed.name.contains(".")}
then
	awk -F '\t' '(\$1~ /^(chr)?[${contig}]\$/)'  ${bed} > TMP/jeter.bed
else
	awk -F '\t' '(\$1~ /^(chr)?[${contig}]\$/) {printf("%s\t0\t%s\\n",\$1,\$2);}' '${fai}' > TMP/jeter.bed
fi

test -s TMP/jeter.bed

bcftools view  -O u --regions-file TMP/jeter.bed --samples-file ${samples}  --trim-unseen-allele --trim-alt-alleles "${vcf}" |\\
	bcftools view -O u -i 'AC>0' |\\
	bcftools annotate -O u -x 'INFO/DP' |\\
	bcftools stats --exons TMP/exons.txt.gz \\
		--fasta-ref '${fasta}' \\
		--samples-file ${samples} > TMP/stats.txt

        
mkdir -p "multiqc.${sex}"
        multiqc  --no-ansi \
                 --title "QC per for ${sex} BED: ${bed.name} contig:${contig}"  \
                 --comment "${vcf.name}"  \
                 --force \
                 --outdir "multiqc.${sex}.${bed.name}.${contig}" \\
		 TMP/stats.txt



zip -9 -r "multiqc.${sex}.${bed.name}.${contig}.zip" "multiqc.${sex}.${bed.name}.${contig}" 
"""
}



process MULTIQC_MALE_VS_FEMALE {
label "process_single"
tag "${directory} ${sample2sex.name} ${contig}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
        tuple val(sex),val(bed_name),val(contig),path(directory),path(sample2sex)
output:
        path("multiqc.${bed_name}.${contig}.M_vs_F.zip"),optional:true,emit:output
script:
"""
hostname 1>&2
mkdir -p TMP/DATA TMP/OUT2

export LC_ALL=en_US.utf8


java -jar \${HOME}/packages/jvarkit/dist/jvarkit.jar multiqcpostproc --sample2collection "${sample2sex}" -o TMP/OUT2 ${directory}/*_data/

find TMP/OUT2 -type f -name "*.json" > TMP/jeter.list

if test -s TMP/jeter.list
then

        mkdir -p "multiqc.M_vs_F"
        multiqc  --filename  "multiqc_report.html" --no-ansi \
                        --title "Males vs females"  \
                        --comment "Males vs Females"  \
                        --force \
                        --outdir "multiqc.${bed_name}.${contig}.M_vs_F" \
                        --file-list TMP/jeter.list

        zip -9 -r multiqc.${bed_name}.${contig}.M_vs_F.zip multiqc.${bed_name}.${contig}.M_vs_F
fi
"""
}
