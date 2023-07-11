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

def gazoduc = gazoduc.Gazoduc.getInstance(params).putGenomeId()


gazoduc.make("bams","NO_FILE").
	description("file containing the path to multiple bam files").
	required().
	existingFile().
	put()


gazoduc.make("beds","NO_FILE").
	description("Path to a list of bed files. if undefined, the REF will be split into parts").
	put()


gazoduc.make("pedigree","NO_FILE").
	description("Optional path to a pedigree").
	put()


gazoduc.make("mapq",10).
	description("mapping quality").
	put()


params.disableFeatures=""


include {runOnComplete;hasFeature;parseBoolean;getKeyValue} from '../../../modules/utils/functions.nf'
include {GATK4_HAPCALLER_GVCFS_01} from '../../../subworkflows/gatk/gatk4.hapcaller.gvcfs.01.nf'
include {COLLECT_TO_FILE_01 as COLLECT2FILE1; COLLECT_TO_FILE_01 as COLLECT2FILE2} from '../../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_CONCAT_01} from '../../../subworkflows/bcftools/bcftools.concat.01.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.02.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {SCATTER_TO_BED} from '../../../subworkflows/picard/picard.scatter2bed.nf' addParams(
		OUTPUT_TYPE:"ACGT",
		MAX_TO_MERGE:"1000"
		)
include {LINUX_SPLIT} from '../../../modules/utils/split.nf' addParams(
	suffix:".bed",
	method:"--lines=1"
	)



if( params.help ) {
    gazoduc.usage().
        name("HC").
        description("Haplotype Caller").
        print();
    exit 0
    }
else
    	{
	gazoduc.validate()
        }


workflow  {
	version_ch = Channel.empty()
   
	/* if no bed was specified, split the genome into part */
        if(file(params.beds).name.equals("NO_FILE")) {
		scatter_ch = SCATTER_TO_BED(params.genomes[params.genomeId].fasta) 
		version_ch = version_ch.mix(scatter_ch.version)

		split_ch = LINUX_SPLIT(scatter_ch.bed)
		version_ch = version_ch.mix(split_ch.version)

		beds_ch = split_ch.output.splitText().map{it.trim()}
	} else {
		beds_ch = Channel.fromPath(params.beds)
	}

	vcfs_ch = GATK4_HAPCALLER_GVCFS_01([:], params.genomeId, file(params.bams), beds_ch, file(params.pedigree))
	version_ch = version_ch.mix(vcfs_ch.version)

	file_list_ch = COLLECT2FILE1([:],vcfs_ch.region_vcf.map{T->T[1]}.collect())
	concat_ch = BCFTOOLS_CONCAT_01([:],file_list_ch.output, file("NO_FILE") )
	version_ch = version_ch.mix(concat_ch.version)

	version2_ch = MERGE_VERSION("Calling gatk", version_ch.collect())

	html_ch = VERSION_TO_HTML(version2_ch.version)

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

runOnComplete(workflow)
