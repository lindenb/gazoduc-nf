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

include {moduleLoad;assertFileExists;isBlank} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {DOWNLOAD_SOMALIER} from '../../modules/somalier/somalier.download.nf'
include {SOMALIER_DOWNLOAD_SITES} from '../../modules/somalier/somalier.download.sites.nf'
//CELUI LA est à changer ==> include {VCF_INTER_PED_01} from '../../modules/bcftools/vcf.inter.pedigree.01.nf'



workflow SOMALIER_VCF_01 {
	take:
		meta
		genomeId
		vcf
		pedigree
	main:
		version_ch = Channel.empty()
		exe_ch = DOWNLOAD_SOMALIER([:])
		version_ch = version_ch.mix(exe_ch.version)

		sites_ch = SOMALIER_DOWNLOAD_SITES([:], genomeId)
		version_ch = version_ch.mix(sites_ch.version)
	
		if(!pedigree.name.equals("NO_FILE")) {
			vcf_samples_ch = VCF_INTER_PED_01(["pedigree_type","other"],vcf,pedigree)	
			version_ch = version_ch.mix(vcf_samples_ch.version)
			ped_ch = vcf_samples_ch.pedigree
			}
		else
			{
			ped_ch = pedigree
			}
	
		somalier_ch = APPLY_SOMALIER([:], genomeId , exe_ch.executable ,sites_ch.vcf ,vcf, ped_ch)
		version_ch = version_ch.mix(somalier_ch.version)

		version_ch = MERGE_VERSION("somalier.vcf", version_ch.collect())
	emit:
		version = version_ch
		zip = somalier_ch.zip
	}


process APPLY_SOMALIER {
tag "${vcf.name}"
afterScript "rm -rf extracted TMP"
input:
	val(meta)
	val(genomeId)
	val(somalier)
	val(sites)
	path(vcf)
	path(pedigree)
output:
	path("${params.prefix?:""}somalier.vcf.zip"),emit:zip
	path("version.xml"),emit:version
script:
	def prefix = params.prefix?:""
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
"""
hostname 1>&2
${moduleLoad("bcftools")}
set -o pipefail
set -x

mkdir TMP

bcftools query -f '%CHROM\t%POS0\t%END\\n' "${sites}" > TMP/jeter.bed

## extract samples
if [ ! -z "${pedigree.name.equals("NO_FILE")?"":"Y"}" ] ; then
	# remove header
	grep -v "^fid\tiid"  "${pedigree}" > TMP/jeter.ped
	test -s TMP/jeter.ped

	cut -f 2 TMP/jeter.ped | sort | uniq > TMP/jeter.samples.txt
	test -s TMP/jeter.samples.txt
fi


if [ ! -z "${vcf.name.endsWith(".list")?"Y":""}" ] && [ -f TMP/jeter.samples.txt ] ; then

	bcftools concat --file-list "${vcf.toRealPath()}" --regions-file TMP/jeter.bed -O u  --allow-overlaps --remove-duplicates |\
		bcftools view --samples-file TMP/jeter.samples.txt -O b -o TMP/jeter.bcf

elif [ ! -z "${vcf.name.endsWith(".list")?"Y":""}" ] ; then

	bcftools concat --file-list "${vcf.toRealPath()}" --regions-file TMP/jeter.bed -O b  --allow-overlaps --remove-duplicates -o TMP/jeter.bcf

elif [ -f TMP/jeter.samples.txt ] ; then

	bcftools view --regions-file TMP/jeter.bed --samples-file TMP/jeter.samples.txt -O b -o TMP/jeter.bcf "${vcf.toRealPath()}"

else

	bcftools view --regions-file TMP/jeter.bed -O b -o TMP/jeter.bcf "${vcf.toRealPath()}"

fi

bcftools index TMP/jeter.bcf


${somalier} extract -d extracted \
	--sites "${sites}" \
	-f "${reference}" TMP/jeter.bcf

mkdir "${prefix}somalier.vcf"

find \${PWD}/extracted -type f -name "*.somalier" | xargs ${somalier} relate --output-prefix=${prefix}somalier.vcf/${prefix}vcf ${pedigree.name.equals("NO_FILE")?"":"-p TMP/jeter.ped"}

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">run somalier on vcf</entry>
        <entry key="vcf">${vcf}</entry>
        <entry key="pedigree">${pedigree}</entry>
</properties>
EOF



# may not exist
touch "${prefix}somalier.vcf/${prefix}vcf.groups.tsv"
zip -9 -r "${prefix}somalier.vcf.zip" "${prefix}somalier.vcf"

"""
}
