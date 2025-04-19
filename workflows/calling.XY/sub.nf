workflow MAKE_STATS {
	take:
		fasta
		fai
		gtf
		samplesheet
		bcf_files
		sex
	main:
		ch1 = EXTRACT_SAMPLES(samplesheet,sex)
		ch2 = BCFTOOLS_STATS(fasta, fai, gtf, bcf_files, sex, ch1.samples)
		if(sex.contains("M") && sex.contains("F")) {
			MULTIQC_MALE_VS_FEMALE(ch1.sample2sex,ch2.directory);
			}
	}


process EXTRACT_SAMPLES {
label "process_quick"
tag "${sex}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	path(samplesheet)
	val(sex)
output:
	path("sample2sex.tsv"),emit:sample2sex
	path("samples.txt"),optional:true,emit:samples
script:
	def sexm = sex.contains("M")?"male":"_ignore"
	def sexf = sex.contains("F")?"female":"_ignore"
"""
mkdir -p TMP

cat << '__EOF__' > TMP/jeter.R
data <- read.table("${samplesheet}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
filtered_data <- data[data\$sex == "${sexm}" | data\$sex == "${sexf}", ]
filtered_data <- filtered_data[, c("sample", "sex")]
write.table(filtered_data, file = "sample2sex.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
__EOF__

R --vanilla --no-save < TMP/jeter.R
cut -f1 sample2sex.tsv | sort | uniq > samples.txt

if test ! -s samples.txt
then
	rm samples.txt
fi
"""
}

process BCFTOOLS_STATS {
label "process_quick"
tag "${sex}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	path(fasta)
	path(fai)
	path(gtf)
	path(bcf_files)
	val(sex)
	path(samples)
output:
	path("multiqc.${sex}.zip"),emit:output
	path("multiqc.${sex}"),emit:directory
script:
	def vcf= bcf_files.find{it.name.endsWith("f")}
"""
mkdir -p TMP
set -o pipefail
export LC_ALL=en_US.utf8

gunzip -c "${gtf}" |\\
	awk -F '\t' '(\$3=="exon" && \$1 ~ /^(chr)?[XY]\$/)' |\\
	cut -f 1,4,5 |\\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\\
	bgzip > TMP/exons.txt.gz


bcftools stats --exons TMP/exons.txt.gz \\
	--fasta-ref '${fasta}' \\
	--regions `cut -f 1 "${fai}" | awk '(\$1~ /^(chr)?[XY]\$/)' | paste -sd,` \\
	--samples-file ${samples} \\
	${vcf} > TMP/stats.txt

        
mkdir -p "multiqc.${sex}"
        multiqc  --no-ansi \
                 --title "QC per for ${sex}"  \
                 --comment "${vcf.name}"  \
                 --force \
                 --outdir "multiqc.${sex}" \\
		 TMP/stats.txt



zip -9 -r "multiqc.${sex}.zip" "multiqc.${sex}" 
"""
}



process MULTIQC_MALE_VS_FEMALE {
label "process_quick"
tag "${directory}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
        path(sample2pop)
        path(directory)
output:
        path("multiqc.M_vs_F.zip"),emit:output
script:
"""
hostname 1>&2
mkdir -p TMP/DATA TMP/OUT2

export LC_ALL=en_US.utf8


java -jar \${HOME}/packages/jvarkit/dist/jvarkit.jar multiqcpostproc --sample2collection "${sample2pop}" -o TMP/OUT2 ${directory}/*_data/

find TMP/OUT2 -type f -name "*.json" > TMP/jeter.list


        mkdir -p "multiqc.M_vs_F"
        multiqc  --filename  "multiqc_report.html" --no-ansi \
                        --title "Males vs females"  \
                        --comment "Males vs Females"  \
                        --force \
                        --outdir "multiqc.M_vs_F" \
                        --file-list TMP/jeter.list

        zip -9 -r multiqc.M_vs_F.zip multiqc.M_vs_F
"""
}
