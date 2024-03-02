include {moduleLoad} from '../../modules/utils/functions.nf'

def TAG="CLINVAR"

workflow ANNOTATE_CLINVAR {
	take:
		genomeId
		vcfs /** tuple vcf,vcf_index */
	main:

	
		if(params.genomes[genomeId].containsKey("clinvar_vcf_url")) {
			source_ch = DOWNLOAD_CLINVAR(genomeId)
			annotate_ch = ANNOTATE(source_ch.vcf, source_ch.index,vcfs)
			out1 = annotate_ch.output
			out2 = annotate_ch.count
			}
		else
			{
			out1 = vcfs
			out2 = Channel.empty()
			}
	emit:
		output = out1
		count = out2
	}


process DOWNLOAD_CLINVAR {
afterScript "rm -rf TMP"
memory "3g"
input:
        val(genomeId)
output:
       	path("${TAG}.db.bcf"),emit:vcf
        path("${TAG}.db.bcf.csi"),emit:index
script:
	def genome = params.genomes[genomeId]
       	def url = genome.clinvar_vcf_url
        def whatis="ClinVar aggregates information about genomic variation and its relationship to human health."

"""
hostname 1>&2
${moduleLoad("bcftools jvarkit")}
set -o pipefail
mkdir -p TMP

wget -O - "${url}" |\
	gunzip -c |\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar vcfsetdict -R "${genome.fasta}"  -n SKIP |\
	bcftools sort -T TMP/tmp -O b -o TMP/clinvar.bcf

bcftools view --header-only TMP/clinvar.bcf | grep "^##INFO" | cut -d '=' -f3 | cut -d ',' -f1 | grep -v '^CLN' | awk '{printf("INFO/%s\tCLN_%s\\n",\$1,\$1);}' > TMP/rename.tsv
bcftools annotate --rename-annots TMP/rename.tsv -O b -o TMP/jeter.bcf TMP/clinvar.bcf
mv TMP/jeter.bcf TMP/clinvar.bcf

bcftools index --force TMP/clinvar.bcf

mv TMP/clinvar.bcf ./${TAG}.db.bcf
mv TMP/clinvar.bcf.csi ./${TAG}.db.bcf.csi
"""
}


process ANNOTATE {
tag "${vcf.name}"
afterScript "rm -rf TMP"
input:
	path(clinvar)
	path(clinvar_idx)
	tuple path(vcf),path(vcf_idx),path(bed)
output:
	tuple path("OUTPUT/${TAG}.bcf"),path("OUTPUT/${TAG}.bcf.csi"),path(bed),emit:output
	path("OUTPUT/count.tsv"),emit:count
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP


bcftools annotate -a "${clinvar}" -c "CLNSIG,CLN_ALLELEID" -O b -o TMP/${TAG}.bcf '${vcf}'
bcftools index --force TMP/${TAG}.bcf
                

###  
bcftools query -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/count.tsv
mv TMP OUTPUT
"""
}
