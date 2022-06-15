include {getKeyValue;getModules} from '../../modules/utils/functions.nf'
include {DOWNLOAD_GTF_01} from '../../modules/gtf/download.gtf.01.nf'
include {COLLECT_TO_FILE} from '../../modules/util/xx.nf'

workflow UORF {
	take:
		meta
		reference
		vcf
		bed
	main:
		version_ch  = Channel.empty()

		gtf_ch = DOWNLOAD_GTF_01(meta, reference )
		version_ch = version_ch.mix(gtf_ch.version)

		vcf2bed_ch = vcf2bed(meta,reference)
		version_ch = version_ch.mix(vcf2bed_ch.version)


		uorf_ch = uorf_gtf(meta, reference, gtf_ch)
		VEP_ANNOT_GTF_01(meta,reference,	
	
		COLLECT_TO_FILE()

		annotate_contig = 
	emit:
	}


process uorf_gtf {
tag "${file(gtf).name}"
input:
	val(meta)
	val(reference)
	val(gtf)
output:
	path("uorf.gtf.gz"),emit:gtf
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
set -o pipefail
module load ${getModules("jvarkit")}

gunzip -c "${gtf}" |\
	java -jar -R "${reference}" |\
	LC_ALL=C sort -T . -t '\t' -k1,1 -k4,4 |\
	bgzip > uorg.gtf.gz

tabix -p gff uorg.gtf.gz

#################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}<entry>
	<entry key="name">extract uORF <entry>
	<entry key="input">${gtf}<entry>
	<entry key="reference">${reference}<entry>
</properties>
EOF
"""
}
