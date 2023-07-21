

workflow BEAGLE01 {
take:
	meta
	genomeId
	rows
main:
	version_ch = Channel.empty()
	
	jar_ch = DOWNLOAD_BEAGLE_JAR([:])
	version_ch = version_ch.mix(jar_ch.version)


	gm_ch = DOWNLOAD_GENETIC_MAPS([:], genomeId)
	version_ch = version_ch.mix(jgm_ch.version)

	rows2_ch = rows.combine( gm_ch.output.splitCsv(sep:'\t',header:false) ).
		filter{T->T[0].contig.equals(T[1][0])}.
		map{T->T[0].plus(gmap:T[1][1])}.
		



	if(params.genomes[genomeId].ucsc_name.equals("hg38")) {
		ch1_ch = DOWNLOAD_REF_HG38_MANIFEST([:])
		version_ch = version_ch.mix(ch1_ch.version)

		ch2_ch = DOWNLOAD_REF_HG38_VCF([:], genomeId, ch1_ch.output.splitCsv(sep:'\t',header:false) )
		version_ch = version_ch.mix(ch2_ch.version)

		rows2_ch = rows2_ch.combine( ch2_ch.output.splitCsv(sep:'\t',header:false) ).
			filter{T->T[0].contig.equals(T[1][0])}.
			map{T->T[0].plus(ref:T[1][1])}.

		}
	else
		{
		ch3_ch =  Channel.empty()
		}


	bgl_ch = APPLY_BEAGLE([:], jar_ch.jar, rows)	
	version_ch = version_ch.mix(blg_ch.version)
	
	out_ch = bgl_ch.output.map{T->T[0].plus(phased_vcf:T[1])}
	
emit:
	version = version_ch
	output = out_ch
}


process DOWNLOAD_BEAGLE_JAR {
tag "${params.beagle_jar_url}"
input:
	val(meta)
output:
	path("beagle.jar"),emit:jar
	path("version.xml"),emit:version
script:
def url = params.beagle_jar_url
"""
hostname 1>&2

wget -O beagle.jar "${url}"

cat << EOF > version.xml
EOF
"""
}


process DOWNLOAD_REF_HG38_MANIFEST {
intput:
	val(meta)
output:
	path("urls.tsv"),emit:output
script:
def base="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"
"""
wget -O - "${base}/phased-manifest_July2021.tsv" |\
	cut -f 1 |\
	grep '.vcf.gz\$' |\
	awk '{printf("${base}\t%s\\n",\$0);}' > urls.tsv
"""
}

pprocess DOWNLOAD_REF_HG38_VCF {
tag "${vcf}"
maxForks 1
input:
	val(meta)
	val(genomeId)
	tuple val(base),val(vcf)
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}
wget -O "${vcf}" "${base}/${vcf}"
wget -O "${vcf}.tbi" "${base}/${vcf}.tbi"

bcftools index -s "${vcf}" | tr "\\n" "\\t" > chrom2vcf.txt
echo "\${PWD}/${vcf}" >> chrom2vcf.txt
"""
}


process DOWNLOAD_GENETIC_MAPS {
tag "${genomeId}"
afterScript "rm -rf TMP OUT/plink.README.txt OUT/README.txt"
input:
	val(meta)
	val(genomeId)
output:
	path("chrom2map.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def genome = params.genome[genomeId]
	def url = genome.beagle_genetic_map
"""
hostname 1>&2
mkdir -p TMP

wget -O TMP/jeter.zip "${url}"

(cd TMP && unzip jeter.zip)

mv -v TMP OUT

# change prefix if needed
find \${PWD}/OUT -type f -name "plink*.map" | while read F
do
	java -jar "\${F}" > TMP/jeter.map
	mv -v TMP/jeter.map "\${F}"

done

find \${PWD}/OUT -type f -name "plink*.map" |\
	awk -F '.' '{printf("%s\t%s\\n",\$(NF-2),\$0);}' |\
	java -jar "\${F}" |\
	sort -T OUT -t '\t' -k1,1V > chrom2map.tsv
cat << EOF > version.xml
EOF
"""
}

process APPLY_BEAGLE {
tag ""
afterScript "rm -rf TMP"
input:
	val(meta)
	path(beagle)
	val(row)
output:
	tuple val(row),path("phased.bcf"),emit:output
script:
	if(!row.containsKey("vcf")) throw new IllegalArgumentException("row.vcf is missing")
"""
hostname 1>&2
mkdir -p TMP

ln -s "${row.vcf}" TMP/jeter.vcf.gz


java -jar ${beagle} ref=${row.ref} map=${row.gmap} gt=TMP/jeter.vcf.gz out=TMP/phased

bcftools view -O b -o phased.bcf TMP/phased.vcf.gz
bcftools index phased.bcf
"""
}
