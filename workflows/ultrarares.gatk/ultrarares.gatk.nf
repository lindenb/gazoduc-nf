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

def gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults().putReference()

gazoduc.make("bams","NO_FILE").
        description("File containing the paths to the BAM/CRAMS files. One path per line").
	required().
	existingFile().
        put()

gazoduc.make("vcf","NO_FILE").
        description("initial vcf file or use --bams").
        put()

gazoduc.make("bed","NO_FILE").
        description("limit to that bed file").
        put()

gazoduc.make("mapq",-1).
        description("mapping quality").
	setInteger().
        put()

gazoduc.make("vcf2interval_distance","25mb").
        description("split VCF per region of 'x' size").
        put()

gazoduc.make("gatkjar","/LAB-DATA/BiRD/users/lindenbaum-p/packages/gatk/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar").
        description("path to gatk jar").
	existingFile().
        put()

include {ULTRA_RARES_ITERATION as ITER_FIRST} from './iteration.part.nf'
include {ULTRA_RARES_ITERATION as ITER_SECOND} from './iteration.part.nf'
include {ULTRA_RARES_ITERATION as ITER_THIRD} from './iteration.part.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {runOnComplete;moduleLoad;getVersionCmd} from '../../modules/utils/functions.nf'
include {SIMPLE_PUBLISH_01} from '../../modules/utils/publish.simple.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'


if( params.help ) {
    gazoduc.usage().
        name("Ultrarares").
        desc("ultrares").
        print();
    exit 0
} else {
   gazoduc.validate();
}



workflow {
	ch1 = ULTRA_RARES_GATK(params, params.reference, params.vcf, file(params.bams), file(params.bed) )
	html = VERSION_TO_HTML(params,ch1.version)
	}

runOnComplete(workflow);


workflow ULTRA_RARES_GATK {
take:
	meta
	reference
	vcf
	bams
	bed
main:
	version_ch = Channel.empty()

	splitbams_ch = SPLIT_BAMS([:], bams)
        version_ch = version_ch.mix(splitbams_ch.version)

	iter1_ch  = ITER_FIRST(meta.plus([split_vcf_method:" --variants-count  100000 "]), reference,vcf, splitbams_ch.output1, bed)
        version_ch = version_ch.mix(iter1_ch.version)

	iter2_ch  = ITER_SECOND(meta.plus([split_vcf_method:" --variants-count  10000 "]), reference, iter1_ch.vcf, splitbams_ch.output2, bed)
        version_ch = version_ch.mix(iter2_ch.version)

	iter3_ch  = ITER_THIRD(meta.plus([split_vcf_method:" --variants-count  1000 "]), reference, iter2_ch.vcf, splitbams_ch.output3, bed)
        version_ch = version_ch.mix(iter3_ch.version)

        version_ch = MERGE_VERSION(meta, "rare-gatk", "rare-gatk",version_ch.collect())

emit:
        vcf = iter1_ch.vcf
        index = iter1_ch.index
        version= version_ch

}

process SPLIT_BAMS {
tag "${bams}"
executor "local"
afterScript "rm jeter.list"
input:
	val(meta)
	path(bams)
output:
	path("list1.bams.list"),emit:output1
	path("list2.bams.list"),emit:output2
	path("list3.bams.list"),emit:output3
	path("version.xml"),emit:version
script:
"""
# sort, biggest bams first
xargs -L1 -a "${bams}" stat  --format="%s %n" | LC_ALL=C sort -t ' ' -k1,1nr | cut -d ' ' -f 2 -> jeter.list

test -s "jeter.list"

head -n 10 "jeter.list" > list1.bams.list

test -s list1.bams.list

tail -n +11 "jeter.list" | head -n 100 > list2.bams.list

tail -n +111 "jeter.list" > list3.bams.list


if test ! -s "list2.bams.list" ; then
	head -n 1 "jeter.list" > list2.bams.list
fi

if test ! -s "list3.bams.list" ; then
	head -n 1 "jeter.list" > list3.bams.list
fi

#####
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">split bams into parts for each iteration</entry>
        <entry key="bams">${bams}</entry>
</properties>
EOF
"""
}
