
workflow VCF_TO_BED {
take:
	vcf_ch
main:
	vcf_ch.branch {
		vcflist : it.name.endsWith(".list")
		vcf: it.name.endsWith(".vcf.gz") || it.name.endsWith(".bcf")
		other: true
		}.set{ ch1 }

	ch1.other.map{throw new IllegalArgumentException("unexpected VCF_TO_BED got $it")}

	ch1.vcflist.
		splitText().
		map{file(it.trim())}.
		mix(ch1.vcf).
		map{[it,file(it.toRealPath() +  (it.name.endsWith(".bcf")?".csi":".tbi"))]}.
		set { ch2 }

	ch3 = VCF2BED(ch2)
	ch4 = MERGE(ch3.output)
emit:
	bed = ch4.output
	chromosomes = ch4.chromosomes
	contigs = ch4.contigs
}

process VCF2BED {
label "process_single"
tag "${vcf.name} ${idx.name}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	tuple path(vcf),path(idx)
output:
	path("${vcf.name}.bed"),emit:output
script:
	if(!(vcf.name.endsWith(".vcf.gz") || vcf.name.endsWith(".bcf")))  IllegalArgumentException("not a valid vcf suffix .vcf.gz .bcf")
"""
set -o pipefail
bcftools index -s "${vcf}" |\\
	awk -F '\t' '{printf("%s\t0\t%s\t${vcf.toRealPath()}\t${idx.toRealPath()}\\n",\$1,\$2);}' > ${vcf.name}.bed
"""
}

process MERGE {
label "process_single"
afterScript "rm -rf TMP"
input:
	path("BED/*")
output:
	path("vcf2bed.bed"),emit:output
	path("chroms.txt"),emit:chromosomes
	path("contigs.bed"),emit:contigs
script:
"""
mkdir -p TMP
set -o pipefail
find BED/ -name "*.bed" > TMP/jeter.list
xargs -a TMP/jeter.list cat |\\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter.bed

cut -f 1 TMP/jeter.bed | uniq |\\
	sort | uniq > TMP/chroms.txt

cut -f1,23 TMP/jeter.bed |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n | uniq > TMP/contigs.bed


mv TMP/contigs.bed ./
mv TMP/jeter.bed vcf2bed.bed
mv TMP/chroms.txt ./
"""
}
