include {readContigFile;hasFeature;parseBoolean;getKeyValue} from '../../../modules/utils/functions.nf'
include {GATK4_HAPCALLER_GVCFS_01} from '../../../subworkflows/gatk/gatk4.hapcaller.gvcfs.01.nf'
include {ANNOTATE} from '../../../subworkflows/annotation/annotation.vcf.01.nf'
include {COLLECT_TO_FILE_01 as COLLECT2FILE1; COLLECT_TO_FILE_01 as COLLECT2FILE2} from '../../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_CONCAT_01} from '../../../subworkflows/bcftools/bcftools.concat.01.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'

params.bam=""
params.beds=""
params.bams=""
params.reference=""
params.references=""
params.dbsnp=""
params.prefix=""
params.publishDir=""
params.disableFeatures=""


workflow  {
	version_ch = Channel.empty()
	vcfs_ch = GATK4_HAPCALLER_GVCFS_01(params,params.reference,file(params.bams),Channel.fromPath(params.beds))
	version_ch = version_ch.mix(vcfs_ch.version)

	if(parseBoolean(getKeyValue(params,"annotate","true"))) {
		bedvcf_ch = COLLECT2FILE1(["header":"interval\tvcf"],vcfs_ch.region_vcf.map{T->T[0]+"\t"+T[1]}.collect())
		version_ch = version_ch.mix(bedvcf_ch.version)


		cfg2=readContigFile(params.annotation_config?:"${workflow.projectDir}/../../../confs/annotation.cfg").plus(params)

		annot_vcf= ANNOTATE(cfg2,params.reference,bedvcf_ch.output.splitCsv(header:true,sep:'\t'))
		version_ch = version_ch.mix(annot_vcf.version)

		file_list_ch = COLLECT2FILE2([:],annot_vcf.bedvcf.splitCsv(header:false,sep:'\t').map{T->T[1]}.collect())
		version_ch = version_ch.mix(file_list_ch.version)
		}
	else
		{
		file_list_ch = COLLECT2FILE1([:],vcfs_ch.region_vcf.map{T->T[1]}.collect())
		}
	concat_ch = BCFTOOLS_CONCAT_01([:],file_list_ch.output)
	version_ch = version_ch.mix(concat_ch.version)


	version2_ch = MERGE_VERSION(params, "Calling gatk", "Calling gatk in gvcf mode", version_ch.collect())



	html_ch = VERSION_TO_HTML(params,version2_ch.version)

	PUBLISH(params,concat_ch.vcf,html_ch.html,version2_ch.version)
	}

process PUBLISH {
executor "local"
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	val(meta)
	val(vcf)
	val(html)
	val(xml)
output:
	path("${params.prefix}genotyped.vcf.gz")
	path("${params.prefix}genotyped.html")
	path("${params.prefix}genotyped.xml")
script:
"""
module load bcftools
bcftools view -O z -o "${params.prefix}genotyped.vcf.gz" "${vcf}"
ln -s "${html}" ./${params.prefix}genotyped.html
ln -s "${xml}" ./${params.prefix}genotyped.xml
"""
}


workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}

