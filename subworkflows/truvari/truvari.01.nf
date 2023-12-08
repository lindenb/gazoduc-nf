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

include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {moduleLoad;getVersionCmd} from '../../modules/utils/functions.nf'
include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {BCFTOOLS_CONCAT_01} from '../bcftools/bcftools.concat.01.nf'


workflow TRUVARI_01 {
     take:
        meta /* meta */
       	genomeId /* genome id */
        vcfs /* file containing the path to VCF files . One per line */
     main:
        version_ch = Channel.empty()

	vcf2bed_ch = VCF_TO_BED([:] ,vcfs)
	version_ch = version_ch.mix( vcf2bed_ch.version)

	perctg_ch = PER_CONTIG([:], genomeId, vcf2bed_ch.chromosomes.splitText().map{it.trim()},vcfs)
	version_ch = version_ch.mix(perctg_ch.version)


	x3_ch = COLLECT_TO_FILE_01([:], perctg_ch.vcf.collect())
	version_ch = version_ch.mix(x3_ch.version)


	concat_ch = BCFTOOLS_CONCAT_01([:], x3_ch.output, file("NO_FILE"))
	version_ch = version_ch.mix(concat_ch.version)
		
		
	version_ch = MERGE_VERSION("Truvari",version_ch.collect())

	emit:
		vcf = concat_ch.vcf
		index = concat_ch.index
		version= version_ch
	}

process PER_CONTIG {
    tag "${contig}"
    afterScript "rm -rf TMP"
    conda "${params.conda}/truvari"
    input:
        val meta
	val genomeId
	val contig
	val vcfs
    output:
	path("truvari.bcf"),     emit: vcf
	path("truvari.bcf.csi"), emit: index
        path "version.xml",emit: version
    script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
	def extraCmd = params.truvari.extra_truvari?:""
	def extraBcftools = params.truvari.extra_view_truvari?:""
    """
	hostname 1>&2
	${moduleLoad("bcftools")}
	mkdir -p TMP

	# bug FORMAT pour dragen
	bcftools merge --force-samples --filter-logic '+' --regions "${contig}" --file-list "${vcfs}" -m none -O u |\
		${extraBcftools.isEmpty()?"":"bcftools view -O u ${extraBcftools} |"} \
		bcftools annotate --force -x 'FORMAT/SR' -O z -o TMP/merged.vcf.gz
	bcftools index --tbi TMP/merged.vcf.gz

	# invoke truvari
	truvari collapse ${extraCmd} --reference "${reference}" -i "TMP/merged.vcf.gz" -c TMP/collapsed.vcf.gz |\
		bcftools +fill-tags -O u -- -t AN,AC,AF |\
		bcftools sort -T TMP -O b -o truvari.bcf

	bcftools index truvari.bcf	

	##########################################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Merge VCFs</entry>
		<entry key="contig">${contig}</entry>
		<entry key="truvari.version">\$( truvari --help 2>&1 | grep ^Truvar) </entry>
		<entry key="truvari.args">${extraCmd}</entry>
		<entry key="versions">${getVersionCmd("bcftools")}</entry>
	</properties>
	EOF
    """
   }
