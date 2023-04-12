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

def gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults().putReference()


gazoduc.make("bams","NO_FILE").
	description("file containing the path to multiple bam files").
	required().
	existingFile().
	put()


gazoduc.make("bed","NO_FILE").
	description("build a pangenome for each segment of this bed").
	put()

gazoduc.make("mapq",10).
	description("mapping quality").
	put()


params.conda=""


include {moduleLoad;runOnComplete;hasFeature;parseBoolean;getKeyValue} from '../../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {SAMTOOLS_SAMPLES_01} from '../../../subworkflows/samtools/samtools.samples.01.nf'

workflow {
	PGGB_01(params, params.reference, file(params.bams), Channel.fromPath(params.bed))
	}


if( params.help ) {
    gazoduc.usage().
        name("pggb").
        description("pggb").
        print();
    exit 0
    }
else
    	{
	gazoduc.validate()
        }



workflow PGGB_01 {
        take:
                meta
                reference
                bams
                bed
        main:
                version_ch = Channel.empty()


		samples_bams_ch = SAMTOOLS_SAMPLES_01([:], reference, file("NO_FILE"), bams)
		version_ch = version_ch.mix(samples_bams_ch.version)

		each_sample_bam = samples_bams_ch.output.
				splitCsv(header:false,sep:'\t').
				map{T->[T[0],T[2]]}

		each_interval = bed.splitCsv(header:false,sep:'\t').
			map{T->T[0]+":"+((T[1] as int)+1)+"-"+T[2]}
		consensus_ch = CALL_SAMTOOLS(meta, reference, each_sample_bam.combine(each_interval))

		APPLY_PGGB(meta, consensus_ch.output.groupTuple())
	}

process CALL_SAMTOOLS {
	tag "${sample} in ${interval} bam:${file(bam).name}"
	afterScript "rm -rf TMP"
	input:
		val(meta)
		val(reference)
		tuple val(sample),val(bam),val(interval)
	output:
		tuple val(interval),path("${sample}.fa"),emit:output
		path("version.xml"),emit:version
	script:
	"""
	${moduleLoad("bcftools samtools")}
	mkdir -p TMP
	bcftools mpileup --regions "${interval}" -Ou -f "${reference}" "${bam}" |\
		bcftools call -mv -Ob -o TMP/calls.bcf

	# normalize indels
	bcftools norm -f "${reference}" TMP/calls.bcf -Ob -o TMP/calls.norm.bcf

	bcftools index TMP/calls.norm.bcf

	# apply variants to create consensus sequence
	samtools faidx "${reference}" "${interval}" | bcftools consensus TMP/calls.norm.bcf > TMP/consensus.fa

	echo ">${sample}.${interval}" > ${sample}.fa

	tail -n +2  TMP/consensus.fa >> "${sample}.fa"


	touch version.xml
	"""
	}

process APPLY_PGGB {
tag "${interval} N=${L.size()}"
conda "${params.conda}/PGGB"
cpus 5
input:
	val(meta)
	tuple val(interval),val(L)
script:
	def otherArgs = "-p 90 -s 100"
"""
hostname 1>&2
${moduleLoad("htslib samtools")}
mkdir -p TMP OUT
cat ${L.join(" ")} > TMP/jeter.fa

bgzip -@ ${task.cpus} TMP/jeter.fa
samtools faidx TMP/jeter.fa.gz


pggb -i TMP/jeter.fa.gz \
     -o OUT \
     -n ${L.size()} \
     ${otherArgs} \
     -t ${task.cpus} 1>&2

touch version.xml

"""
}
