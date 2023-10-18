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


include {moduleLoad;runOnComplete} from '../../../modules/utils/functions.nf'
include {MOSDEPTH_DOWNLOAD_01} from '../../../modules/mosdepth/mosdepth.downoad.01.nf'
include {MOSDEPTH_RUN_01} from '../../../modules/mosdepth/mosdepth.run.01.nf'
include {SAMTOOLS_SAMPLES02} from '../../../subworkflows/samtools/samtools.samples.02.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.02.nf'
include {MULTIQC_01} from '../../../modules/multiqc/multiqc.01.nf' 
include {COLLECT_TO_FILE_01} from '../../../modules/utils/collect2file.01.nf'
include {SCATTER_TO_BED} from '../../../subworkflows/picard/picard.scatter2bed.02.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {SOMALIER_BAMS_02} from '../../../subworkflows/somalier/somalier.bams.01.nf'


if( params.help ) {
    exit 0
}


workflow {
	ch1 = BAMSTATS([:],file(params.bams), file(params.bed) , file(params.sample2pop) )
	html = VERSION_TO_HTML(ch1.version)
	}


workflow BAMSTATS {
take:
	meta
	bams
	bed
	sample2pop
main:
	version_ch = Channel.empty()
	toqc_ch = Channel.empty()

        samples_ch = SAMTOOLS_SAMPLES02( [with_header:true,allow_multiple_references: true,allow_duplicate_samples : true], "", bams)
        version_ch = version_ch.mix(samples_ch.version)
	
	rows_ch = samples_ch.output.splitCsv(header:true,sep:"\t")


	rows_ch = rows_ch.
		map{T->T.plus("sample":T.new_sample)}.
		map{T->{
			def genomeId = params.genomes.grep{it.value.fasta.equals(T.reference)}.collect{it.key}.join("");
			if(genomeId.isEmpty()) throw new IllegalStateException("cannot find genomeId for ${T}");
			return T.plus([genomeId:genomeId]);
		}}



	if(params.with_somalier) {

		genomeId_ch = rows_ch.
        	    map{T->T.genomeId}.
		    first()

		som_ch = SOMALIER_BAMS_02(
			[:],
			genomeId_ch,
			rows_ch.combine(genomeId_ch).
				filter{T->T[0].genomeId.equals(T[1])}
				.map{T->T[0]},
			file("NO_FILE") /* pedigree */
			)
                version_ch = version_ch.mix(som_ch.version)
		toqc_ch =  toqc_ch.mix(som_ch.qc)
		}


	if(bed.name.equals("NO_FILE")) {

		/** unique ref */
		refs_ch = rows_ch.
			map{T->[reference:T.reference]}.
			unique()
		acgt_ch = SCATTER_TO_BED(["OUTPUT_TYPE":"ACGT","MAX_TO_MERGE":"1"], refs_ch )
                version_ch = version_ch.mix(acgt_ch.version)

		rows_ch = rows_ch.
			combine(acgt_ch.output).
			filter{T[0].reference.equals(T[1].reference)}.
			map{T->T[0].plus("bed":T[1].scatter_bed)}

		}
	else
		{
		rows_ch = rows_ch.map{T->T.plus("bed":bed)}
		}


	if(params.with_samtools_idxstats) {
		fstats_ch = APPLY_SAMTOOLS_IDXSTAT([:],rows_ch)
		version_ch = version_ch.mix(fstats_ch.version)
		toqc_ch =  toqc_ch.mix(fstats_ch.output.map{T->T[1]})
		}
	

	if(params.with_samtools_flagstat) {
		fstats_ch = APPLY_SAMTOOLS_FLAGSTAT([:],rows_ch)
		version_ch = version_ch.mix(fstats_ch.version)
		toqc_ch =  toqc_ch.mix(fstats_ch.output.map{T->T[1]})
		}

	if(params.with_samtools_stats) {
		ststats_ch = APPLY_SAMTOOLS_STATS([:],rows_ch)
		version_ch = version_ch.mix(ststats_ch.version)

		if(sample2pop.name.equals("NO_FILE")) {
			toqc_ch =  toqc_ch.mix(ststats_ch.output.map{T->T[1]})
			}
		else
			{
			statsperpop_ch = SAMTOOLS_STATS_PER_POP([:], ststats_ch.output.map{T->T[0].sample+"\t"+T[1]}.collect(), sample2pop)
			version_ch = version_ch.mix(statsperpop_ch.version)

			toqc_ch =  toqc_ch.mix(ststats_ch.map{T->T[1]})
			}


		}
 
	if(params.with_CollectWgsMetrics) {
		ch1_ch = APPLY_COLLECT_WGS_METRICS([:],rows_ch)
		version_ch = version_ch.mix(ch1_ch.version)
		toqc_ch =  toqc_ch.mix(ch1_ch.output.map{T->T[1]})
		}

	if(params.with_mosdepth) {
                mosdepth_ch = MOSDEPTH_DOWNLOAD_01([:])
                version_ch = version_ch.mix(mosdepth_ch.version)
		
		
                ch2 = MOSDEPTH_RUN_01([:], mosdepth_ch.executable,rows_ch)
                version_ch = version_ch.mix(ch2.version)
		

		 toqc_ch =  toqc_ch.mix(
			ch2.output.splitCsv(sep:'\t',header:true).
				flatMap{T->[T.globaldist,T.regiondist,T.summary]}.
				filter{T->!(T==null || T.equals("NO_FILE") || T.equals(".") || T.equals(""))}
			)
		/* 
                merge_ch = MERGE_MOSDEPTH_SUMMARY(meta,ch2.summary.map{T->T[0].sample+"\t"+T[0].bam+"\t"+T[0].reference+"\t"+T[1]}.collect())
                version_ch = version_ch.mix(merge_ch.version)
		*/
		}

         multiqc_ch = MULTIQC_01([extra:" --fullnames "], toqc_ch.collect())
         version_ch = version_ch.mix(multiqc_ch.version)


         version_ch = MERGE_VERSION("BAM stats",version_ch.collect())
emit:
	version = version_ch
}


runOnComplete(workflow);


process APPLY_SAMTOOLS_IDXSTAT {
tag "${row.sample}"
cpus 1
input:
      	val(meta)
        val(row)
output:
       	tuple val(row),path("*.idxstats.txt"),emit:output
        path("version.xml"),emit:version
when:
	row.bam.endsWith(".bam") // idx stats is just slow with CRAM...
script:
	if(!row.containsKey("genomeId")) throw new IllegalStateException("cannot find genomeId for ${row}");
	def genomeId = row.genomeId
        def prefix = row.prefix?:(params.prefix?:"") + genomeId + "."
"""
hostname 1>&2
${moduleLoad("samtools")}

samtools view idxstats "${row.bam}" > "${row.sample}.${prefix}idxstats.txt"

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">samtools idxstats</entry>
        <entry key="sample">${row.sample}</entry>
        <entry key="bam">${row.bam}</entry>
        <entry key="samtools.version">\$(samtools version | head -n1 )</entry>
</properties>
EOF
"""
}



process APPLY_SAMTOOLS_FLAGSTAT {
tag "${row.sample}"
afterScript "rm -rf TMP"
cpus 1
input:
      	val(meta)
        val(row)
output:
       	tuple val(row),path("*.flags.txt"),emit:output
        path("version.xml"),emit:version
script:
	if(!row.containsKey("genomeId")) throw new IllegalStateException("cannot find genomeId for ${row}");
	def genomeId = row.genomeId
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
        def prefix = row.prefix?:(params.prefix?:"") + genomeId + "."
	def extraST= (row.bed.name.equals("NO_FILE")?"":" -M --regions-file \"${row.bed}\" ")
"""
hostname 1>&2
${moduleLoad("samtools")}
set -o pipefail
mkdir -p TMP

samtools view ${extraST} -O BAM --uncompressed --reference "${reference}" "${row.bam}" |\
samtools flagstats '-' > TMP/jeter.flags.txt


mv -v "TMP/jeter.flags.txt" "${row.sample}.${prefix}flags.txt" 	


##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">samtools flagstats</entry>
        <entry key="sample">${row.sample}</entry>
        <entry key="bam">${row.bam}</entry>
        <entry key="samtools.version">\$(samtools version | head -n1 )</entry>
</properties>
EOF
"""
}


process APPLY_SAMTOOLS_STATS {
tag "${row.sample}"
afterScript "rm -rf TMP"
cpus 1
input:
      	val(meta)
        val(row)
output:
       	tuple val(row),path("*.stats.txt"),emit:output
        path("version.xml"),emit:version
script:
	if(!row.containsKey("genomeId")) throw new IllegalStateException("cannot find genomeId for ${row}");
	def genomeId = row.genomeId
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
        def prefix = row.prefix?:(params.prefix?:"") + genomeId + "."
        def extra = params.samtools.stats.args?:"--remove-dups "
	def extraST= (row.bed.name.equals("NO_FILE")?"":" -M --regions-file \"${row.bed}\" ")
"""
hostname 1>&2
${moduleLoad("samtools")}
set -o pipefail
mkdir -p TMP

samtools view ${extraST} -O BAM --uncompressed --reference "${reference}" "${row.bam}" |\
samtools stats ${extra}  --ref-seq "${reference}" --reference "${reference}" > TMP/jeter.stats.txt

awk '/^# The command line/ {printf("# sample : ${row.sample}\\n");} {print}' TMP/jeter.stats.txt > "TMP/jeter2.stats.txt"
mv -v TMP/jeter2.stats.txt TMP/jeter.stats.txt


## FIX ME , BUG IN MULTIQC https://github.com/ewels/MultiQC/issues/1971

awk -F '\t' 'BEGIN {F=0;} {OFS="\t"; if(\$1=="SN") {if(\$2=="filtered sequences:") {F=int(\$3);\$3="0";} if(\$2=="reads unmapped:" && F>0) {\$3=int(\$3)+F;} } print;}' TMP/jeter.stats.txt > "TMP/jeter2.stats.txt"
mv -v TMP/jeter2.stats.txt TMP/jeter.stats.txt



mv -v "TMP/jeter.stats.txt" "${row.sample}.${prefix}stats.txt" 	


##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">samtools stats</entry>
        <entry key="sample">${row.sample}</entry>
        <entry key="bam">${row.bam}</entry>
        <entry key="samtools.version">\$(samtools version | head -n1 )</entry>
</properties>
EOF
"""
}



process APPLY_COLLECT_WGS_METRICS {
tag "${row.sample}"
afterScript "rm -rf TMP"
memory "3g"
input:
        val(meta)
        val(row)
output:
        tuple val(row),path("*.collect_wgs_metrics.txt"),emit:output
        path("version.xml"),emit:version
script:
       	if(!row.containsKey("genomeId")) throw new IllegalStateException("cannot find genomeId for ${row}");
        def genomeId = row.genomeId
        def genome = params.genomes[genomeId]
        def reference = genome.fasta
        def prefix = row.prefix?:(params.prefix?:"") + genomeId + "."
"""
hostname 1>&2
${moduleLoad("gatk/0.0.0")}
set -o pipefail
mkdir -p TMP

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"  CollectWgsMetrics \
	I=${row.bam} \
	O=TMP/jeter.txt \
	R=${reference} \
	${params.gatk.collectWgsMetrics.args}

mv -v TMP/jeter.txt "${row.sample}.${prefix}.collect_wgs_metrics.txt"
 

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">gatk CollectWgsMetrics</entry>
	<entry key="extra-params">${params.gatk.collectWgsMetrics.args}</entry>
</properties>
EOF
"""
}

