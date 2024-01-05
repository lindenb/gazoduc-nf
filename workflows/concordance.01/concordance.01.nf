/*

Copyright (c) 2024 Pierre Lindenbaum

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

def gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults().putReference()


gazoduc.make("vcf","NO_FILE").
	description("input vcf").
	existingFile().
	required().
	put()

gazoduc.make("bed","NO_FILE").
	description("restrict to that bed").
	put()



include {runOnComplete} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {SAMPLES_IN_VCF_01} from '../../modules/bcftools/samples.in.vcf.01.nf'
include {moduleLoad;isBlank;isHg38;isHg19} from '../../modules/utils/functions.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'

if( params.help ) {
    gazoduc.usage().
	name("concordance01").
	desc("Concordance").
	print();
    exit 0
} else {
   gazoduc.validate();
}



workflow {
	ch = CONCORDANCE_01(params, params.reference, params.vcf , file(params.bed))
	//html = VERSION_TO_HTML(params,ch.version)
	//SIMPLE_PUBLISH_01(params, Channel.empty().mix(html.html).mix(ch.version).mix(ch.zip).collect())
	}

runOnComplete(workflow);

workflow CONCORDANCE_01 {
    take:
	    meta
	    reference
	    vcf
	    bed
    main:
		version_ch = Channel.empty()
		
		sn_ch = SAMPLES_IN_VCF_01([:],vcf)
		version_ch = version_ch.mix(sn_ch.version)
	
		rgn_ch = MAKE_INTERVALS([:],reference, bed)		
		version_ch = version_ch.mix(rgn_ch.version)

		ch_sns_1 = sn_ch.samples.splitText().map{it.trim()}
		ch_sns_2 = ch1_sns

		each_pair_ch = ch_sns_1.combine(ch_sns_2).
			filter{T->!T[0].equals(T[1])}
		
		conc_ch = CALC_CONCORDANCE([:],vcf, rgn_ch.intervals, each_pair_ch)
		version_ch = version_ch.mix(conc_ch.version)

		//version_ch = MERGE_VERSION(meta, "Concordance", "Concordance", version_ch.collect())
    emit:
	 //   version = version_ch
    }

process MAKE_INTERVALS {
tag ${bed.name}"
memory "2g"
input:
	val(meta)
	val(reference)
	path(bed)
output:
	path("intervals.interval_list"),emit:intervals
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("gatk/0.0.0")}
mkdir -p TMP

if ${bed.name.equals("NO_FILE")} ; then

	awk -F '\t' '{printf("%s\t0\t%d\\n",\$1,\$2);}' '${reference}.fai' > TMP/jeter.bed

else

	cp "${bed}" TMP/jeter.bed

fi

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" 
		BedToIntervalList \
        	SEQUENCE_DICTIONARY=${reference} \
		O=TMP/jeter.interval_list UNIQUE=true I=TMP/jeter.bed

mv TMP/jeter.interval_list intervals.interval_list

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
	<entry key="versions">${getVersionCmd("awk gatk")}</entry>
</properties>
EOF
"""
}

process CALC_CONCORDANCE {
tag "${sn1} ${sn2}"
memory "2g"
input:
	val(meta)
	val(vcf)
	path(intervals)
	tuple val(sn1),val(sn2)
output:
	tuple val(sn1),
		val(sn2),
		path("${sn1}_vs_${sn2}.genotype_concordance_summary_metrics"),
		path("${sn1}_vs_${sn2}.genotype_concordance_detail_metrics"),
		path("${sn1}_vs_${sn2}.genotype_concordance_contingency_metrics")
	path("version.xml"),emit:version
script:
	def min_dp=20
	def min_gp=30
"""
hostname 1>&2
${moduleLoad("gatk/0.0.0")}
mkdir -p TMP

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" GenotypeConcordance \
        CALL_VCF=${vcf} \
        CALL_SAMPLE=${sn1} \
        TRUTH_VCF=${vcf} \
        TRUTH_SAMPLE=${sn2} \
        O=${sn1}_vs_${sn2} \
        INTERVALS=${intervals} \
	MIN_DP=${min_dp} \
	MIN_GQ=${min_gq}

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
	<entry key="versions">${getVersionCmd("gatk")}</entry>
</properties>
EOF
"""
}
