
params.genomeId= "";
params.vcf="NO_FILE"
params.samples="NO_FILE"
params.bed = "NO_FILE"

include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf'
include {BEAGLE01} from '../../subworkflows/beagle/beagle.01.nf'
include {moduleLoad;runOnComplete;getVersionCmd} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_CONCAT_01} from '../../subworkflows/bcftools/bcftools.concat.01.nf'



workflow {
	DO_BEAGLE([:], params.genomeId, file(params.vcf),file(params.bed),file(params.samples) )
	}

runOnComplete(workflow)


workflow DO_BEAGLE {
	take:
		meta
		genomeId
		vcf
		bed
		samples
	main:
		version_ch= Channel.empty()

		vcf2bed_ch = VCF_TO_BED([with_header:false], vcf)
		version_ch = version_ch.mix(vcf2bed_ch.version)

		if(bed.name.equals("NO_FILE")) {
			bed_ch = vcf2bed_ch.bed 
			}
		else
			{
			inter_ch = INTERSECT_BED([:], vcf2bed_ch.bed, bed)	
			version_ch = version_ch.mix(inter_ch.version)

			bed_ch = inter_ch.bed
			}


		intervals2_ch = bed_ch.splitCsv(sep:'\t',header:false).map{T->[
			contig:T[0],
			interval:T[0]+":"+((T[1] as int) +1)+"-"+T[2],
			vcf:T[3],
			samples: samples
			]} 


		bgl_ch = BEAGLE01([:], genomeId, intervals2_ch)
		version_ch = version_ch.mix(bgl_ch.version)


		x3_ch = COLLECT_TO_FILE_01([:], bgl_ch.output.map{T->T.phased_vcf}.collect())
		version_ch = version_ch.mix(x3_ch.version)

		x4_ch = BCFTOOLS_CONCAT_01([:],x3_ch.output,file("NO_FILE"))
		version_ch = version_ch.mix(x4_ch.version)


                version_ch = MERGE_VERSION("beagle",version_ch.collect())
	emit:
		version = version_ch
		vcf = x4_ch.vcf
	}


process INTERSECT_BED {
tag "${vcfbed.name} / ${userbed.name}"
afterScript "rm -rf TMP"
input:
        val(meta)
        path(vcfbed)
        path(userbed)
output:
        path("intersect.bed"),emit:bed
        path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bedtools")}
set -o pipefail

mkdir -p TMP

sort -t '\t' -T TMP -k1,1 -k2,2n "${vcfbed}" | cut -f1-4 >  TMP/a.bed
sort -t '\t' -T TMP -k1,1 -k2,2n "${userbed}" | cut -f1,2,3 | bedtools merge >  TMP/b.bed

# 4th column contains the path to the VCF
bedtools intersect -a TMP/a.bed -b TMP/b.bed |\
        sort -t '\t' -T TMP -k1,1 -k2,2n | uniq > intersect.bed

#######################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">intersection between vcf bed and user bed</entry>
        <entry key="version">${getVersionCmd("bedtools")}</entry>
</properties>
EOF
"""
}


