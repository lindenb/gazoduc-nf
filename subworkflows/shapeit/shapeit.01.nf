include {MERGE_VERSION} from '../../modules/version/version.merge.nf'

/*


*/

workflow SHAPEIT_01 {
take:
	meta
	genomeId
main:
	version_ch = Channel.empty()

	shapeit_ch = DOWNLOAD_SHAPEIT([:],genomeId)
	version_ch = version_ch.mix(shapeit_ch.version)

	pc_ch = PHASE_COMMON([:],shapeit_ch.phase_common, row)

	LIGATE_COMMON([:], pc_ch.ligate , pc_ch.map{T->[ T[0].substring(0,T[0].indexOf(':')), T[1] ]}.groupTuple() )

	version_ch = MERGE_VERSION("ShapeIt",version_ch.collect())
emit:
	version_ch = version_ch

}

workflow SHAPEIT_COMMON {
	take:
		meta
		genomeId
		vcf
	emit:
		version_ch = version_ch
		
	}


workflow QC_PER_CHUNK {
	take:
		meta
		genomeId
		vcf
	main:
		


}


process PHASE_COMMON {
cpus 10
input:
	val(meta)
	val(phase_common)
	val(row)
output:
	tuple val("${row.interval}"),path("phased.bcf"),emit:output
script:
	def maf=row.maf?:0.001
"""
hostname 1>&2

mkdir -p TMP
${phase_common} --input '${row.vcf}' \
	--map '${row.map}' \
	--output TMP/jeter.bcf \
	--thread '${task.cpu}' \
	--log TMP/jeter.log \
	--filter-maf ${maf} \
	--region '${row.interval}'

bcftools index -f  --threads '${task.cpu}' TMP/jeter.bcf

mv TMP/jeter.bcf ./phased.bcf
mv TMP/jeter.bcf.csi ./phased.bcf.csi
"""
}

process LIGATE_COMMON  {
input:
	val(meta)
	val(ligate)
	tuple (L)
script:
"""
hostname 1>&2
mkdir -p TMP

cat << EOF > TMP/jeter.list
${L.join("\n")}
EOF

${ligate} --input TMP/jeter.list --output jeter.bcf --thread ${task.cpus} --index
"""
}


JOBID1=$(dx run app-swiss-army-knife -iin="/Bulk/Whole\ genome\ sequences/Whole\ genome\ GraphTyper\ joint\ call\ pVCF/ukb23352_c${CHR}_b${CHUNK}_v1.vcf.gz" 
-iin="/Bulk/Whole\ genome\ sequences/Whole\ genome\ GraphTyper\ joint\ call\ pVCF/ukb23352_c${CHR}_b${CHUNK}_v1.vcf.gz.tbi" 
-iin="${ODIR0}/pass_qc_chr${CHR}.bcf" 
-iin="${ODIR0}/pass_qc_chr${CHR}.bcf.csi" 
-icmd="


	JOBID0=$(dx run app-swiss-army-knife -icmd="

tabix /mnt/project/Bulk/Whole\ genome\ sequences/Whole\ genome\ GraphTyper\ joint\ call\ pVCF/QC/qc_metrics_graphtyper_v2.7.1_qc.tab.gz chr${CHR} | awk 'BEGIN{print \"##fileformat=VCFv4.2\n##contig=<ID=chr${CHR}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\"}(10^$14<1 && 10^$14>0){AF=10^\$14;MAF=(AF>0.5?1-AF:AF);ExcHet=\$17/(2*(\$16+\$17+\$18)*(AF)*(1-AF));}{if (\$7>0.5 && \$15 <15012 && ExcHet >=0.5 && ExcHet <= 1.5){ split(\$3, array, \":\"); print \$1\"\t\"\$2\"\t\"\$4\"\t\"array[3]\"\t\"array[4]\"\t.\t.\t.\t.\"}}' | bcftools view -Ob -o pass_qc_chr${CHR}.bcf && bcftools index --threads 2 -f pass_qc_chr${CHR}.bcf" --tag filter_tab --tag tabix --tag swiss_army_knife --instance-type mem1_ssd1_v2_x2 --folder="${ODIR0}" --name qc_step03a --priority normal -y | tail -n1 | cut -d" " -f3)

bcftools annotate --regions-file ${bed} -x '^INFO/AC,^INFO/AN,^FORMAT/GT' -O u '${vcf}' |\
bcftools view -f PASS -Ou |\
bcftools norm -m -any -Ou |\
bcftools view -i 'ALT!="*"' -Ob -o genotyped.bcf

bcftools index -f genotyped.bcf

bcftools isec -c none -n=2 -w1 --regions-file '${bed}' 'genotyped.bcf' 'pass_qc_chr${CHR}.bcf' -Ob -o pass_qc_ukb23352_c${CHR}_${START}_${END}_v1.bcf && bcftools index -f pass_qc_ukb23352_c${CHR}_${START}_${END}_v1.bcf && bcftools +fill-tags -r chr${CHR}:${START}-${END} -Ou pass_qc_ukb23352_c${CHR}_${START}_${END}_v1.bcf -- -t HWE -S /mnt/project/data/ukb_wgs/unphased/qc/support/UKB_samples_with_WGS_british_single_gbr.txt | bcftools view -G -e \"INFO/HWE_GBR < 1e-30\" -Ob -o pass_hwe_chr${CHR}_${START}_${END}.bcf && bcftools index -f pass_hwe_chr${CHR}_${START}_${END}.bcf && bcftools isec -c none -n=2 -w1 -r chr${CHR}:${START}-${END} pass_qc_ukb23352_c${CHR}_${START}_${END}_v1.bcf pass_hwe_chr${CHR}_${START}_${END}.bcf -Ob -o ukb23352_c${CHR}_${START}_${END}_v1.bcf && bcftools index -f ukb23352_c${CHR}_${START}_${END}_v1.bcf && rm -f pass* split_*" --tag qc --tag filter_tab --tag filter_hwe --instance-type mem1_ssd1_v2_x2 --folder="${ODIR1}" --depends-on ${JOBID0} --name qc_chr${CHR}_chunk${CHUNK} --priority normal -y | tail -n1 | cut -d" " -f3)
