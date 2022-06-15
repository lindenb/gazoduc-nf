process SNPEFF_ANNOT_01 {
afterScript "rm -f jeter.bcf"
input:
	val(meta)
	val(vcf)
	val(snpEffCfg)
	val(dbName)
output:
	path("snpeff.bcf"),emit:vcf
	path("snpeff.bcf.csi"),emit:csi
script:
"""
hostname 1>&2
set -o pipefail
module load ${getModules("snpEff bcftools")}

java -Djava.io.tmpdir=. -jar "$${SNPEFF_JAR}" eff \
	-config "${snpEffCfg}" -interval "jeter.bed" \
	-nodownload -noNextProt -noMotif -noInteraction -noLog -noStats -chr chr -i vcf -o vcf ${dbName} "${vcf}" |\
	bcftools view -o jeter.bcf -O b

bcftools index jeter.bcf

mv jeter.bcf "snpeff.bcf"
mv jeter.bcf.csi "snpeff.bcf.csi"
"""
}
