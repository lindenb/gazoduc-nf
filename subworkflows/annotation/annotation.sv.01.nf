/*

Copyright (c) 2023 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
The MIT License (MIT)
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

def gazoduc = gazoduc.Gazoduc.getInstance()

gazoduc.make("sv_annot_gtf","NO_FILE").
        description("Location of TABIX indexed GTF for annotation").
	setBoolean().
        put()


include {getVersionCmd;jvarkit;isHg19;isHg38;isBlank;hasFeature;moduleLoad;getGnomadGenomePath;getGnomadExomePath} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {DOWNLOAD_GNOMAD_SV_01} from '../gnomad/download_gnomad_sv.01.nf'
include {DOWNLOAD_DGV_01} from '../../modules/dgv/download.dgv.01.nf'
include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf'
include {JENSENLAB_DISEASES_01} from '../jensenlab/diseases.01.nf'
include {BCFTOOLS_CONCAT_01} from '../bcftools/bcftools.concat.01.nf'
include {CONCAT_FILES_01} from '../../modules/utils/concat.files.nf'


workflow ANNOTATE_SV_VCF_01 {
	take:
		meta
		reference
		vcf
	main:
		version_ch = Channel.empty()
	
		gtf_file = file(params.sv_annot_gtf)


		gnomad_ch = DOWNLOAD_GNOMAD_SV_01(meta, reference)
		version_ch = version_ch.mix(gnomad_ch.version)		

		dgv_ch = DOWNLOAD_DGV_01([:], reference) 
		version_ch = version_ch.mix(dgv_ch.version)

		ensemblreg_ch = DOWNLOAD_ENSEMBL_REG([:], reference)
		version_ch = version_ch.mix(ensemblreg_ch.version)

		diseases_ch = JENSENLAB_DISEASES_01([:], reference, gtf_file )

		ch1_ch  = VCF_TO_BED(meta,vcf)
		version_ch = version_ch.mix(ch1_ch.version)

		ch2_ch = ch1_ch.bed.splitCsv(header:false,sep:'\t').
			map{T->[T[0],T[3]]}


		db2_ch = CONCAT_FILES_01(meta, Channel.empty().mix(diseases_ch.output).collect())
		version_ch = version_ch.mix(db2_ch.version)
	

		ann_ch= ANNOTATE(meta, reference, gnomad_ch.bed, dgv_ch.bed, ensemblreg_ch.output, db2_ch.output, gtf_file, ch2_ch)
		version_ch = version_ch.mix(ann_ch.version)
	
		ch5_ch = COLLECT_TO_FILE_01(meta, ann_ch.vcf.collect())
		version_ch = version_ch.mix(ch5_ch.version)

		ch6_ch = BCFTOOLS_CONCAT_01(meta,ch5_ch.output)
		version_ch = version_ch.mix(ch6_ch.version)


		version_ch = MERGE_VERSION(meta, "svannotation", "VCF SV annotation", version_ch.collect())
	emit:
		version= version_ch
		vcf = ch6_ch.vcf
		index = ch6_ch.index
	}


process DOWNLOAD_ENSEMBL_REG {
input:
	val(meta)
	val(reference)
output:
	path("ensembl.reg.gff.gz"),emit:output
	path("version.xml"),emit:version
	path("ensembl.reg.gff.gz.tbi"),emit:index
script:
	def hg19url = "https://ftp.ensembl.org/pub/grch37/current/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.gff.gz"
	def hg38url = "https://ftp.ensembl.org/pub/release-109/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz"
	def url=isHg19(reference)? hg19url :(isHg38(reference)? hg38url:"")
"""
hostname 1>&2
${moduleLoad("htslib")}

if ${url.isEmpty()} ; then

	touch ensembl.reg.gff

else

	wget -O - "${url}" | gunzip -c | LC_ALL=C sort -T . -t '\t' -k1,1 -k4,4n > ensembl.reg.gff
	test -s ensembl.reg.gff

fi

bgzip ensembl.reg.gff
tabix -f -p gff ensembl.reg.gff.gz

cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">Download ensembl ref</entry>
	<entry key="url"><a>${url}</a></entry>
	<entry key="version">${getVersionCmd("wget bgzip tabix")}</entry>
</properties>
EOF
"""
}

process ANNOTATE {
tag "${file(vcf).name} ${contig}"
memory "3g"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
	val(gnomad_sv_bed)
	val(dgv_bed)
	val(ensembl_reg)
	path(annotations_files)
	val(gtf)
	tuple val(contig),val(vcf)
output:
	path("output.bcf"),emit:vcf
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("jvarkit bedtools bcftools")}
mkdir -p TMP

bcftools view "${vcf}" "${contig}" |\
	java  -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar \
			vcfsvannotator \
			--dgv '${dgv_bed}' \
			--gnomad '${gnomad_sv_bed}' \
			--ensemblreg '${ensembl_reg}' \
			--gtf '${gtf}' |\
	bcftools view -O u -o TMP/output.bcf

	sed -e 's/,GENE_NAME,DOID,/,-,-,/' "${annotations_files}" | while IFS='\t' read DATABASE COLS HEADER
        do
                test -f "\${DATABASE}"
                test -f "\${HEADER}"
                test -f "\${DATABASE}.tbi"

                bcftools annotate --force --annotations "\${DATABASE}" \
                        -h "\${HEADER}" \
                        -c "\${COLS}" `echo "\${COLS}" | tr "," "\\n" | awk '(\$1!="." && \$1!="-" && NR>3) {S="unique"; printf(" --merge-logic %s:%s ",\$1,S);}'  ` -O u -o TMP/jeter2.bcf TMP/output.bcf
                mv TMP/jeter2.bcf TMP/output.bcf
        done


bcftools index TMP/output.bcf

mv TMP/output.bcf ./
mv TMP/output.bcf.csi ./

cat <<- EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">Annotate contig</entry>
	<entry key="vcf">${vcf}</entry>
	<entry key="contig">${contig}</entry>
	<entry key="version">${getVersionCmd("bcftools")}</entry>
</properties>
EOF
"""
}
