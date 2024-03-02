include {moduleLoad} from '../../modules/utils/functions.nf'

def TAG="SNPEFF"

workflow ANNOTATE_SNPEFF {
	take:
		genomeId
		vcfs /** tuple vcf,vcf_index */
	main:
		source_ch =  BUILD_SNPEFF(genomeId)
		annotate_ch = ANNOTATE(genomeId,source_ch.output, vcfs)
	emit:
		output = annotate_ch.output
		count = annotate_ch.count
}


/** build snpeff Database from gtf */
process BUILD_SNPEFF {

afterScript "rm -f org.tsv genes.tsv snpeff.errors"
memory "10g"
input:
        val(genomeId)
output:
       	path("snpEff.config"),emit:output
script:
        def genome = params.genomes[genomeId]
        def dbName = genomeId
        def reference = genome.fasta
        def gtf = genome.gtf
"""
hostname 1>&2
${moduleLoad("snpEff")}
set -o pipefail

mkdir -p "data/${dbName}"
ln -s "${reference}" "data/${dbName}/sequences.fa"

cp "${gtf}"  data/${dbName}/genes.gtf.gz
gunzip data/${dbName}/genes.gtf.gz

# write snpEff contig
cat << EOF > snpEff.config
data.dir = \${PWD}/data/
${dbName}.genome = Human
${dbName}.reference = ${reference}
EOF

# build database
java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar \${SNPEFF_JAR}  build -gtf22 -v "${dbName}" 2> snpeff.errors

rm snpeff.errors

test -s "data/${dbName}/snpEffectPredictor.bin"

rm  data/${dbName}/genes.gtf
"""
}

process ANNOTATE {
tag "${vcf.name}"
afterScript "rm -rf TMP"
memory '3g'
input:
	val(genomeId)
	path(config)
	tuple path(vcf),path(vcf_idx),path(bed)
output:
	tuple path("OUTPUT/${TAG}.bcf"),path("OUTPUT/${TAG}.bcf.csi"),path(bed),emit:output
	path("OUTPUT/count.tsv"),emit:count
script:
"""
hostname 1>&2
${moduleLoad("snpEff bcftools")}
mkdir -p TMP

bcftools view '${vcf}' -O v |\
java  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar "\${SNPEFF_JAR}" eff -config '${config}' \
                                -nodownload -noNextProt -noMotif -noInteraction -noLog -noStats -chr chr -i vcf -o vcf '${genomeId}' > TMP/jeter1.vcf

bcftools sort --max-mem '${task.memory.giga}G' -T TMP/tmp -O b -o TMP/${TAG}.bcf TMP/jeter1.vcf
bcftools index TMP/${TAG}.bcf

bcftools query -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/count.tsv
mv TMP OUTPUT
"""
}
