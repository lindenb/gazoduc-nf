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
include {SAMTOOLS_SAMPLES02} from '../../subworkflows/samtools/samtools.samples.02.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {isHg19;isHg38;moduleLoad} from '../../modules/utils/functions.nf'
include {BCFTOOLS_CONCAT_01} from '../../subworkflows/bcftools/bcftools.concat.01.nf'


workflow BCFTOOLS_CALL_01 {
	take:
		meta
		genomeId
		bams
		beds
	main:
		version_ch = Channel.empty()

		sn_ch = SAMTOOLS_SAMPLES02(["with_header":false,"allow_multiple_references":true], genomeId, bams)
		version_ch = version_ch.mix(sn_ch.version)

		each_bed = beds.splitText().map{T->T.trim()}

		ref2bams = sn_ch.output.splitCsv(header:false,sep:'\t').
			map{T->[T[3],T[2]]}.
			groupTuple()

		split_per_ref_ch = SPLIT_PER_REF([:], ref2bams)
		version_ch = version_ch.mix(split_per_ref_ch.version)

		ref_bams_bed = split_per_ref_ch.output.splitCsv(header:false,sep:'\t').combine(each_bed)


		call_ch = CALL_INTERVAL([:], genomeId, ref_bams_bed)
		version_ch = version_ch.mix(call_ch.version)

		merge_ch = BCFTOOL_MERGE_BED([:], call_ch.output.groupTuple())
		version_ch = version_ch.mix(merge_ch.version)

		file_list_ch = COLLECT_TO_FILE_01([:],merge_ch.output.collect())
		version_ch = version_ch.mix(file_list_ch.version)

		concat_ch = BCFTOOLS_CONCAT_01([:],file_list_ch.output)
		version_ch = version_ch.mix(concat_ch.version)

                version_ch = MERGE_VERSION("Calling samtools", version_ch.collect())
			
	emit:
		version = version_ch
		vcf = concat_ch.vcf
	}


process SPLIT_PER_REF {
executor "local"
tag "${genomeId} N=${L.size()}"
input:
	val(meta)
	tuple val(genomeId),val(L)
output:
	path("output.txt"),emit:output
	path("version.xml"),emit:version
script:
	def n = params.num_bams_per_call
"""
hostname 1>&2
set -o pipefail
mkdir OUT

cat << EOF | split --lines=${n} --additional-suffix=.list - OUT/cluster.
${L.join("\n")}
EOF

find \${PWD}/OUT -type f -name "cluster.*.list" | awk '{printf("${genomeId}\t%s\\n",\$0);}'  > output.txt

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">split bam list</entry>
        <entry key="genomeId">${genomeId}</entry>
        <entry key="n">${n}</entry>
</properties>
EOF
"""
}

process CALL_INTERVAL {
tag "${genomeId2} ${bams} ${bed}"
afterScript  "rm -rf TMP"
input:
	val(meta)
	val(genomeId)
	tuple val(genomeId2),val(bams),val(bed)
output:
	tuple val(bed),path("output.bcf"),emit:output
	path("version.xml"),emit:version
script:
	def mapq = meta.mapq?:20
	def t= task.cpus?"--threads ${task.cpus}":""
	def genome = params.genomes[genomeId2]
	def reference2 = genome.fasta
	def ploidy = isHg19(reference2)?"--ploidy GRCh37":(isHg38(reference2)?"--ploidy GRCh38":"")
"""
hostname 1>&2
${moduleLoad("bcftools jvarkit")}
set -o pipefail
mkdir -p TMP


if [ "${reference}" != "${reference2}" ] ; then
	## TODO
	false
else
	cp "${bed}" TMP/jeter.bed
fi



bcftools mpileup --redo-BAQ -a 'FORMAT/AD' -a 'FORMAT/DP' -a 'INFO/AD' \
		--fasta-ref "${reference2}" ${t} \
		--regions-file "TMP/jeter.bed" -q ${mapq} -O u --bam-list ${bams} -O u -o TMP/jeter.bcf

bcftools call  --keep-alts ${ploidy} ${t} \
		--multiallelic-caller --variants-only --output-type b -o "TMP/jeter2.bcf" TMP/jeter.bcf

mv TMP/jeter2.bcf TMP/jeter.bcf

bcftools index ${t} TMP/jeter.bcf

if [ "${genomeId}" != "${genomeId2}" ] ; then
	## TODO
	false	
fi

mv TMP/jeter.bcf output.bcf
mv TMP/jeter.bcf.csi output.bcf.csi


#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">call bcftools mpileup and call</entry>
        <entry key="src.reference">${reference2}</entry>
        <entry key="dest.reference">${reference}</entry>
	<entry key="bcftools.version">\$(bcftools version | head -n 2 | paste -s)</entry>
</properties>
EOF
"""
}

process BCFTOOL_MERGE_BED {
tag "${bed} ${L.size()}"
afterScript "rm -rf TMP"
input:
	val(meta)
	tuple val(bed),val(L)
output:
	path("output.bcf"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}
set -o pipefail
mkdir TMP

cat << EOF > TMP/jeter.list
${L.join("\n")}
EOF


#
# use the first sample of each VCF to be sure that samples will be ordered the same way 
# for the downstream bcftools concat
#
cat TMP/jeter.list | while read V
do
	bcftools query -l "\$V" | head -n 1 | tr "\n" "," >> TMP/jeter2.list
	echo "\${V}" >> TMP/jeter2.list
done

sort -t, -k1,1 TMP/jeter2.list | cut -d, -f2 > TMP/jeter3.list


bcftools merge --file-list TMP/jeter3.list -O b -o TMP/jeter.bcf
bcftools index TMP/jeter.bcf

mv TMP/jeter.bcf output.bcf
mv TMP/jeter.bcf.csi output.bcf.csi

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">merge vcf for one bed</entry>
        <entry key="bed">${bed}</entry>
	<entry key="bcftools.version">\$(bcftools version | head -n 2 | paste -s)</entry>
</properties>
EOF
"""
}
