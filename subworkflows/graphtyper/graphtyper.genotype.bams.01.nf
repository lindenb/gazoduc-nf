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



include {GRAPHTYPER_DOWNLOAD_01} from '../../modules/graphtyper/graphtyper.download.01.nf'
include {GRAPHTYPER_GENOTYPE_01} from '../../modules/graphtyper/graphtyper.genotype.01.nf'
//include {SCATTER_TO_BED} from '../../subworkflows/picard/picard.scatter2bed.02.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {isBlank;moduleLoad;getVersionCmd} from '../../modules/utils/functions.nf'
include {MOSDEPTH_DOWNLOAD_01} from '../../modules/mosdepth/mosdepth.downoad.01.nf'
include {COLLECT_TO_FILE_01 as COLLECT_TO_FILE_X1} from '../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_CONCAT_01} from '../../subworkflows/bcftools/bcftools.concat.01.nf'
include {MOSDEPTH_RUN_01} from '../../modules/mosdepth/mosdepth.run.01.nf'
//include {SAMTOOLS_DEPTH_01} from '../../subworkflows/samtools/samtools.depth.01.nf'
include {SAMTOOLS_SAMPLES02} from '../../subworkflows/samtools/samtools.samples.02.nf'
include {CONCAT_FILES_01} from '../../modules/utils/concat.files.nf'


workflow GRAPHTYPER_GENOTYPE_BAMS_01 {
	take:
		meta
		genomeId
		bams
		bed
		sample2depth
	main:
		version_ch = Channel.empty()
		
		concat_bed_ch = MERGE_BED([:],genomeId, bed)
		version_ch = version_ch.mix(concat_bed_ch.version)

		windows_bed_ch = MAKE_WINDOWS_50KB([:], genomeId, concat_bed_ch.output)
		version_ch = version_ch.mix(windows_bed_ch.version)
		

		bams_ch = SAMTOOLS_SAMPLES02(["with_header":false,"allow_multiple_references":false,"allow_duplicate_samples":false],genomeId,bams)

		unknowndp_ch = EXTRACT_SAMPLES_WITH_UNKNOWN_DEPTH([:],sample2depth,bams_ch.output)
		version_ch = version_ch.mix(unknowndp_ch.version)

		
		summary_ch = Channel.empty()
		if(!sample2depth.name.equals("EMPTY")) {
			summary_ch = summary_ch.mix(unknowndp_ch.knowndp)
			}

		if( params.depth_method.equals("mosdepth") ) {

			mosdepth_input_ch  = unknowndp_ch.bams4depth.
				splitCsv(sep:'\t',header:false).
				map{T->[
					"sample":T[0],
					"bam":T[1],
					"genomeId":genomeId,
					"mapq": params.mapq
					]}.
				combine(concat_bed_ch.output).
				map{T->T[0].plus("bed":T[1])}

			mosdepth_ex = MOSDEPTH_DOWNLOAD_01([:])

			mosdepth_ch = MOSDEPTH_RUN_01([:], mosdepth_ex.executable, mosdepth_input_ch)

			merge_mosdepth_ch = MERGE_MOSDEPTH([:], mosdepth_ch.output.map{T->T[1]}.collect() )
			version_ch = version_ch.mix(merge_mosdepth_ch.version)

			summary_ch = summary_ch.mix(merge_mosdepth_ch.output)
			}
		/*
		else if((params.depth_method.equals("samtoolsdepth")) {
			stdp_ch = SAMTOOLS_DEPTH_01([:],genomeId,file("NO_FILE"),bams, bed0)
			version_ch = version_ch.mix(stdp_ch.version)
			
			convert_ch = CONVERT_TO_MOSDEPTH_SUMMARY([:],stdp_ch.output)

			summary= convert_ch.output
			} */
		else
			{
			throw new IllegalArgumentException("unknown params.depth_method")
			}
		
		x1_ch = COVERAGE_DIVIDE_READLENGTH([:], genomeId, summary_ch.collect())
		version_ch = version_ch.mix(x1_ch.version)
		

		executable_ch = GRAPHTYPER_DOWNLOAD_01([:])
		version_ch = version_ch.mix(executable_ch.version)

		
		each_bed_ch = windows_bed_ch.output.splitText().map{it.trim()}
	
		x2_ch = GRAPHTYPER_GENOTYPE_01([:], executable_ch.executable,x1_ch.output.combine(each_bed_ch).map{T->[
			"bams":T[0],
			"avg_cov_by_readlen":T[1],
			"bed":T[2],
			"genomeId":genomeId
			]})

		version_ch = version_ch.mix(x2_ch.version)
		
		x3_ch = CONCAT_FILES_01(["suffix":".list","concat_n_files":50,"downstream_cmd":""],x2_ch.output.map{T->T[1]}.collect())
		version_ch = version_ch.mix(x3_ch.version)


		x4_ch = BCFTOOLS_CONCAT_01([:],x3_ch.output,file("NO_FILE"))
		version_ch = version_ch.mix(x4_ch.version)

		version_ch = MERGE_VERSION("genotype bams with graphtyper", version_ch.collect())
	emit:
		version = version_ch
		vcf = x4_ch.vcf
	}


process MERGE_BED {
executor "local"
tag "${bed.name}"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(genomeId)
	path(bed)
output:
	path("concat.bed"),emit:output
	path("version.xml"),emit:version
script:
"""
set -o pipefail
mkdir -p TMP 
${moduleLoad("bedtools jvarkit")}
cut -f1,2,3 "${bed}"|\
	sort -T . -t '\t' -k1,1 -k2,2n |\
	bedtools merge > TMP/concat.bed

mv TMP/concat.bed ./

test -s concat.bed


cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">concatenate bed file+merge.</entry>
	<entry key="version">${getVersionCmd("bedtools")}</entry>
</properties>
EOF
"""
}


process MAKE_WINDOWS_50KB {
tag "${bed.name}"
memory "3g"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(genomeId)
	path(bed)
output:
	path("windows.list"),emit:output
	path("version.xml"),emit:version
script:
	def w=50000
	def s=w-10
	def reference = params.genomes[genomeId].fasta
"""
set -o pipefail
mkdir -p TMP BEDS
${moduleLoad("bedtools jvarkit")}

bedtools makewindows -b '${bed}' -w "${w}" -s "${s}"|\
	awk -F '\t' '(int(\$2)<int(\$3))' |\
	sort -T . -t '\t' -k1,1V -k2,2n > TMP/windows.bed


java  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar bedcluster \
		-R "${reference}" \
		--size '${params.bedcluster_size}' \
		-o BEDS TMP/windows.bed

find \${PWD}/BEDS -type f -name "*.bed"  > windows.list

test -s concat.bed
test -s windows.list

sleep 10

cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">concatenate bed file+merge. Make windows of 50Kb</entry>
	<entry key="version">${getVersionCmd("bedtools")}</entry>
</properties>
EOF
"""
}


process EXTRACT_SAMPLES_WITH_UNKNOWN_DEPTH {
tag "${sample2depth} ${bams}"
executor "local"
input:
	val(meta)
	path(sample2depth)
	path(bams)
output:
	path("knowndepth.tsv"),emit:knowndp
	path("bams4depths.tsv"),emit:bams4depth
	path("version.xml"),emit:version
script:

if(sample2depth.name.equals("NO_FILE") )

"""
touch knowndepth.tsv

cut -f2,3 '${bams}' > bams4depths.tsv

cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">join know/unknown depth</entry>
</properties>
EOF
"""

else

"""
join -t '\t' -1 1 -2 1 -o '2.2,2.3,1.2' \
	<(sort -T . -t '\t' -k1,1 '${sample2depth}' | grep -v '^#' | uniq) \
	<(sort -T . -t '\t' -k1,1 '${bams}') > knowndepth.tsv


join -t '\t' -1 1 -2 1 -v 2 -o '2.2,2.3' \
	<(sort -T . -t '\t' -k1,1 '${sample2depth}' | uniq) \
	<(sort -T . -t '\t' -k1,1 '${bams}') > bams4depths.tsv

cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">join know/unknown depth</entry>
</properties>
EOF
"""
}



process MERGE_MOSDEPTH {
tag "N=${L.size()}"
executor "local"
input:
	val(meta)
	val(L)
output:
	path("summary.tsv"),emit:output
	path("version.xml"),emit:version
"""
mkdir -p TMP

cat << EOF > TMP/jeter1.txt
${L.join("\n")}
EOF

touch TMP/jeter2.txt

xargs -a TMP/jeter1.txt cat | grep -v '^sample' | while read SN	REF	BAM	_globaldist	_regiondist	SUMMARY	_perbase	_regions
do
        echo -n "\${SN}\t\${BAM}" >> TMP/jeter2.txt
        awk -F '\t' 'BEGIN{C="";} (\$1=="total") {if(C=="") {C=\$4;}} (\$1=="total_region") {C=\$4;} END{printf("\t%s",C);}' "\${SUMMARY}" >> TMP/jeter2.txt
        echo >> TMP/jeter2.txt
done

mv TMP/jeter2.txt "summary.tsv"

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">merge mosdepth summary</entry>
        <entry key="count">${L.size()}</entry>
</properties>
EOF
"""
}

process CONVERT_TO_MOSDEPTH_SUMMARY {
executor "local"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(L)
output:
	path("summary.tsv"),emit:output
	path("version.xml"),emit:version
script:
"""
echo -e "sample\tbam\tref\tcov\tcovr" > jeter.tsv
awk -F '\t' '/^#/ {next;} {printf("%s\t%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,\$6,\$6);}' '${depths}' >> jeter.tsv

sleep 10
mv jeter.tsv summary.tsv
"""
}

process COVERAGE_DIVIDE_READLENGTH {
tag "N=${L.size()}"
input:
	val(meta)
	val(genomeId)
	val(L)
output:
	tuple path("bams.txt"),path("avg_cov_by_readlen.txt"),emit:output
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
	def n_reads = meta.n_reads?:1000
"""
hostname 1>&2
${moduleLoad("samtools")}

cat ${L.join(" ")} | while read SN BAM COV
do
	samtools view -F 3844 -T "${reference}" "\${BAM}" | head -n "${n_reads}" |\
		awk -F '\t' -vCOV=\${COV} 'BEGIN{T=0.0;N=0;} {N++;T+=length(\$10)} END{print COV/(N==0?100:T/N);}' >> avg_cov_by_readlen.txt
	echo "\${BAM}" >> bams.txt
done

cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">create list of bams and list of coverage/read-length</entry>
</properties>
EOF
"""
}
