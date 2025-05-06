

workflow  {
	contigs_ch = Channel.of("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")

	ch1 = DOWNLOAD_GNOMAD(contigs_ch)
	
	ch1bis = MERGE_GNOMAD(ch1.output.collect())
	
	cadd_names=Channel.of("whole_genome_SNVs")
	
	ch2 = DOWNLOAD_CADD(contigs_ch.combine(cadd_names))
	
	ch3 = MERGE_CADD(ch2.output.groupTuple())

	ch5 = DOWNLOAD_VARIANT_CATALOG()
	ch6 = PAR_BED()

	MERGE_CONFIGS(ch3.output.
		mix(ch5.output).
		mix(ch6.output).
		mix(ch1bis.output).
		collect())
	}


process DOWNLOAD_GNOMAD {
tag "${contig}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
label "process_quick"
input:
	val(contig)
output:
	path("gnomad.${contig}.txt"),emit:output
script:
	def pop = params.gnomad_population
"""
mkdir -p TMP
set -o pipefail

wget -O - "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr${contig}.vcf.bgz" | \\
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%${pop}\\n' |\\
	gzip > TMP/gnomad.${contig}.tsv.gz

mkdir -p "${params.outdir}/TMP"

mv TMP/gnomad.${contig}.tsv.gz "${params.outdir}/TMP"
echo "${params.outdir}/TMP/gnomad.${contig}.tsv.gz" > gnomad.${contig}.txt
"""
}

process MERGE_GNOMAD {
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
label "process_quick"
input:
	path("TSV/*")
output:
	path("gnomad.flag"),emit:output
script:
"""
set -o pipefail
mkdir -p TMP
cat TSV/*.txt | LC_ALL=C sort -V > TMP/jeter.list
xargs -a TMP/jeter.list gunzip -c | bgzip > TMP/jeter.tsv.gz

tabix -s 1 -b 2 -e 2 TMP/jeter.tsv.gz
mkdir -p "${params.outdir}/GNOMAD"
mv TMP/jeter.tsv.gz "${params.outdir}/GNOMAD/gnomad.genomes.v4.1.sites.tsv.gz"
mv TMP/jeter.tsv.gz.tbi "${params.outdir}/GNOMAD/gnomad.genomes.v4.1.sites.tsv.gz.tbi"
touch gnomad.flag
rm -v ${params.outdir}/TMP/gnomad.*.tsv.gz
"""
}


process DOWNLOAD_CADD {
afterScript "rm -rf TMP"
tag "${name} ${contig}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
label "process_quick"
input:
	tuple val(contig),val(name)
output:
	tuple val(name),path("${name}.${contig}.txt"),emit:output
script:
"""
mkdir -p TMP
set -o pipefail
tabix "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/${name}.tsv.gz" "${contig}" |\
	awk -F '\t' '/^#/{${contig.equals("1")?"print;":""}next;} {printf("chr%s\\n",\$0);}' |\\
	gzip > TMP/jeter.tsv.gz

rm -vf *.tbi

mkdir -p "${params.outdir}/TMP"
mv TMP/jeter.tsv.gz "${params.outdir}/TMP/tmp.cadd.${name}.${contig}.tsv.gz"
echo "${params.outdir}/TMP/tmp.cadd.${name}.${contig}.tsv.gz" > "${name}.${contig}.txt"
"""
}

process MERGE_CADD {
tag "${name}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
label "process_medium"
input:
        tuple val(name),path("TSV/*")
output:
        path("cadd.${name}.flag"),emit:output
script:
"""
mkdir -p TMP
set -o pipefail
cat TSV/*.txt | LC_ALL=C sort -V > TMP/jeter.list
xargs -a TMP/jeter.list gunzip -c  | bgzip > TMP/jeter.tsv.gz

tabix -s 1 -b 2 -e 2 TMP/jeter.tsv.gz

mkdir -p "${params.outdir}/CADD/GRCh38"
mv TMP/jeter.tsv.gz "${params.outdir}/CADD/GRCh38/${name}.tsv.gz"
mv TMP/jeter.tsv.gz.tbi "${params.outdir}/CADD/GRCh38/${name}.tsv.gz.tbi"


touch "cadd.${name}.flag"
rm -fv ${params.outdir}/TMP/tmp.cadd.${name}.*
"""
}


process DOWNLOAD_VARIANT_CATALOG {
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
label "process_quick"
output:
	path("xhunter.flag"),emit:output
script:
"""
mkdir -p TMP   
wget -O TMP/jeter.json "https://raw.githubusercontent.com/Clinical-Genomics/reference-files/refs/heads/master/rare-disease/disease_loci/ExpansionHunter-v5.0.0/variant_catalog_grch38.json"

mkdir -p "${params.outdir}/XHUNTER"
mv TMP/jeter.json "${params.outdir}/XHUNTER/variant_catalog_grch38.json"

touch "xhunter.flag"
"""
}

process PAR_BED {
executor "local"
output:
	path("par.flag"),emit:output
script:
"""
mkdir -p "${params.outdir}/PAR"

cat << EOF | tr ":-" "\t" > ${params.outdir}/PAR/par.bed
chrX:10001-2781479
chrX:155701383-156030895
chrY:10001-2781479
chrY:56887903-57217415
EOF

touch par.flag
"""
}

process MERGE_CONFIGS {
executor "local"
input:
	path("CFGS/*")
output:
	path("resources.config"),emit:output
script:
"""
cat << EOF > resources.config
params {
	gnomad_af = "${params.outdir}/GNOMAD/gnomad.genomes.v4.1.sites.tsv.gz"
	gnomad_af_idx = "${params.outdir}/GNOMAD/gnomad.genomes.v4.1.sites.tsv.gz.tbi"
	known_dbsnp = "/LAB-DATA/GLiCID/projects/BiRD_resources/species/human/cng.fr/hs38me/hs38me_all_chr_dbsnp-146.vcf.gz"
	known_dbsnp_tbi = "/LAB-DATA/GLiCID/projects/BiRD_resources/species/human/cng.fr/hs38me/hs38me_all_chr_dbsnp-146.vcf.gz.tbi"
	variant_catalog = "${params.outdir}/XHUNTER/variant_catalog_grch38.json"
	par_bed = "${params.outdir}/PAR/par.bed"
	vep_cache = "/LAB-DATA/GLiCID/projects/BiRD_resources/apps/vep"
}
EOF
"""
}
