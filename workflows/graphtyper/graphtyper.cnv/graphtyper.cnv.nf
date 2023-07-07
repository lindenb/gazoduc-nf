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
nextflow.enable.dsl=2

def gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults()

gazoduc.build("vcf","NO_FILE").
	desc("file containing the SV variants to genotype").
	existingFile().
	required().
	put()

gazoduc.make("bams","NO_FILE").
        description("File containing the paths to the BAM/CRAMS files. One path per line").
        required().
        existingFile().
        put()

gazoduc.make("bed","NO_FILE").
        description("Retstict to that bed").
        put()


gazoduc.make("extraBcfTools","").
        description("extra arguments to add to bcftools where splitting vcf (could be option -e 'blblabala' )").
        put()


gazoduc.make("split_vcf_method","--vcf-count 1000").
        description("split initial VCF. Argument for jvarkit/splitnvariants.").
        put()


include {GRAPHTYPER_DOWNLOAD_01} from '../../../modules/graphtyper/graphtyper.download.01.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.02.nf'
include {COLLECT_TO_FILE_01} from '../../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_CONCAT_01} from '../../../subworkflows/bcftools/bcftools.concat.01.nf'
include {JVARKIT_VCF_SPLIT_N_VARIANTS_01} from '../../../subworkflows/jvarkit/jvarkit.vcfsplitnvariants.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {runOnComplete;moduleLoad;getVersionCmd} from '../../../modules/utils/functions.nf'

if( params.help ) {
    gazoduc.usage().
	name("GraphTyper CNV").
	desc("Genotype CNV/SV using graphtyper").
	print();
    exit 0
} else {
   gazoduc.validate();
}



workflow {
	c1_ch = GRAPHTYPER_CNV_01(file(params.bams),file(params.vcf), file(params.bed))
	html = VERSION_TO_HTML(c1_ch.version)	
	}

runOnComplete(workflow);


process PUBLISH {
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	path(html)
	path(version)
output:
	path(zip)
when:
	!params.getOrDefault("publishDir","").trim().isEmpty()
script:
"""
echo "publishing ${zip}"
"""
}




workflow GRAPHTYPER_CNV_01 {
	take:
		bams
		vcf
		bed
	main:
		version_ch = Channel.empty()
				
		executable_ch = GRAPHTYPER_DOWNLOAD_01()
		version_ch = version_ch.mix(executable_ch.version)

		splitvcf_ch = JVARKIT_VCF_SPLIT_N_VARIANTS_01(vcf, bed)
		version_ch = version_ch.mix(splitvcf_ch.version)

		each_vcf = splitvcf_ch.output.splitText().map{it.trim()}

		x2_ch = GENOTYPE_CNV(
			executable_ch.executable, 
			bams,
			each_vcf
			)
		
		x3_ch = COLLECT_TO_FILE_01([:],x2_ch.output.map{T->T[1]}.collect())
		version_ch = version_ch.mix(x3_ch.version)

		x4_ch = BCFTOOLS_CONCAT_01([:],x3_ch.output,file("NO_FILE"))
		version_ch = version_ch.mix(x4_ch.version)

		version_ch = MERGE_VERSION("graptyper",version_ch.collect())
	emit:
		version = version_ch
		vcf = x4_ch.vcf
	}


process GENOTYPE_CNV {
tag "${file(vcf).name}"
cpus 5
afterScript "rm -rf results TMP sv_results"
errorStrategy 'retry'
maxRetries 3
input:
	val(graphtyper)
	path(bams)
	val(vcf)
output:
	path("genotyped.bcf"),emit:output
	path("version.xml"),emit:version
script:
	def reference = params.genomes[params.genomeId].fasta
"""
hostname 1>&2
${moduleLoad("bcftools bedtools")}
mkdir -p TMP

export TMPDIR=\${PWD}/TMP

bcftools query -f '%CHROM\t%POS0\t%END\\n' '${vcf}' | sort -T TMP -t '\t' -k1,1 -k2,2n | bedtools merge |\
	awk -F '\t' '{printf("%s:%d-%s\\n",\$1,int(\$2)+1,\$3);}' > TMP/jeter.intervals

# prevent empty intervals
if test ! -s TMP/jeter.intervals ; then
	awk '(NR==1){printf("%s:1-2\\n",\$1);}' "${reference}.fai" > TMP/jeter.intervals
fi


${graphtyper} genotype_sv \
	"${reference}" \
	"${vcf}" \
	--force_no_copy_reference \
	--force_use_input_ref_for_cram_reading \
	--region_file TMP/jeter.intervals \
	--sams=${bams} \
	--threads=${task.cpus}

find \${PWD}/sv_results/ -type f -name "*.vcf.gz" | grep -v '/input_sites/' > TMP/vcf.list

bcftools concat --file-list TMP/vcf.list \
	--allow-overlaps --remove-duplicates \
	--threads ${task.cpus} -O u |\
	bcftools sort -T TMP -O b -o "genotyped.bcf"

bcftools index --threads ${task.cpus} "genotyped.bcf"


##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">run graphtyper</entry>
        <entry key="bams">${bams}</entry>
        <entry key="reference">${reference}</entry>
        <entry key="version">${getVersionCmd("bcftools bedtools")}</entry>
</properties>
EOF
"""
}

