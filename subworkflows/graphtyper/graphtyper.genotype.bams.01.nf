/*

Copyright (c) 2022 Pierre Lindenbaum

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
include {SCATTER_TO_BED} from '../../subworkflows/picard/picard.scatter2bed.nf'
include {SAMTOOLS_FASTQ_01} from '../../modules/samtools/samtools.collate.fastq.01.nf'
include {MAP_BWA_01} from '../../subworkflows/mapping/map.bwa.01.nf'
include {MERGE_VERSION as MERGE_VERSIONA; MERGE_VERSION as MERGE_VERSIONB} from '../../modules/version/version.merge.nf'
include {isBlank;moduleLoad} from '../../modules/utils/functions.nf'


workflow GRAPHTYPER_GENOTYPE_BAMS_01 {
	take:
		meta
		reference
		bams
		bed
	main:
		version_ch = Channel.empty()
		
		each_interval_ch = bed.splitCsv(header:false,sep:'\t').
			map{T->[
			"bams":bams,
			"reference":reference,
			"interval":T[0]+":"+((T[1] as int)+1)+"-"+T[2]]}.
			]}
	
	emit:
		version = version_ch
		bams = remap_ch.bams
	}

workflow REMAP_ONE {
take:
	meta
	reference_in
	reference_out
	bed
	sample_bam
main:
	version_ch = Channel.empty()

	bamr_ch = BAM_AND_REGIONS(meta, reference_in, bed, sample_bam )
	version_ch = version_ch.mix(bamr_ch.version)

	unmapped_ch = sample_bam.map{T->[
		"sample": T[0],
		"bam": T[1],
		"reference": reference_in,
		"interval": "*"
		]}

	rows_ch = bamr_ch.output.splitCsv(header:false,sep:'\t').map{T->[
		"sample": T[0],
		"bam": T[1],
		"reference": reference_in,
		"bed": T[2]
		]}

	fastq_ch = SAMTOOLS_FASTQ_01(meta,rows_ch.mix(unmapped_ch))
	version_ch = version_ch.mix(fastq_ch.version)

	map1_ch = fastq_ch.output.splitCsv(header:true,sep:'\t')

	umap1_ch = map1_ch.flatMap{T->[
		[[T.sample,"R1"],T.unpairedR1],
		[[T.sample,"R2"],T.unpairedR2]			
		]}.groupTuple()
	unmap2_ch = SORT_UNPAIRED_FASTQ(meta, umap1_ch)
	version_ch = version_ch.mix(unmap2_ch.version)

	single_ch = map1_ch.map{T->[
		"sample" : T.sample,
		"R1" : T.other,
		"reference" : reference_out
		]}.filter{T->!isEmptyGz(T.R1)}
	
	unmap3_R1_ch = unmap2_ch.output.filter{it[1].equals("R1")}.map{T->[T[0],T[2]]}
	unmap3_R2_ch = unmap2_ch.output.filter{it[1].equals("R2")}.map{T->[T[0],T[2]]}
	unmap4_ch = JOIN_UNPAIRED(meta,unmap3_R1_ch.join(unmap3_R2_ch))
	version_ch = version_ch.mix(unmap4_ch.version)

	r1r2_ch= map1_ch.map{T->[
		"sample":T.sample,
		"R1":T.R1,
		"R2":T.R2
		]}.filter{T->!(isEmptyGz(T.R1) || isEmptyGz(T.R2))}
	
	bam2_ch = unmap4_ch.output.map{T->[
                "sample":T[0],  
                "R1":T[1],
		"interleaved":true
                ]}.filter{T->!isEmptyGz(T.R1)}
	
	single2_ch = unmap4_ch.single.map{T->[
                "sample":T[0],
                "R1":T[1]
                ]}.filter{T->!isEmptyGz(T.R1)}
	
	bam1_ch = MAP_BWA_01(meta,reference_out,r1r2_ch.
		mix(bam2_ch).
		mix(single_ch).
		mix(single2_ch))
	version_ch = version_ch.mix(bam1_ch.version)

	version_ch = MERGE_VERSIONB(meta, "Remap", "Remap bam on another reference", version_ch.collect())
emit:
	version= version_ch.version
	bams = bam1_ch.bams
}


process BAM_AND_REGIONS {
tag "${sample}"
afterScript "rm -f jeter.intervals contigs.tsv"
input:
	val(meta)
	val(reference)
	val(bed)
	tuple val(sample),val(bam)
output:
	path("output.tsv"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("samtools jvarkit")}
touch contigs.tsv

# find contigs with at least one read

awk '{printf("%s:%d-%s\\n",\$1,int(\$2)+1,\$3);} END{printf("*\\n");}' "${bed}"  > jeter.intervals


cat jeter.intervals | while read R
do
	samtools view --reference "${reference}" "${bam}" "\${R}" | head -n1 |\
		awk -vR=\"\${R}\" '{printf("%s\\n",R);}' >> contigs.tsv
done

mkdir BED
awk -F '[:\\-]' '{printf("%s\t%d\t%s\\n",\$1,int(\$2)-1,int(\$3));}' "contigs.tsv" |\
	java -jar \${JVARKIT_DIST}/bedcluster.jar --reference "${reference}" -o BED --size "${meta.collateSize?:"100mb"}" 

find \${PWD}/BED -type f -name "*.bed" |\
	awk '{printf("${sample}\t${bam}\t%s\\n",\$0);}' > output.tsv


##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">find regions in BAM where there is a leat one read</entry>
        <entry key="sample">${sample}</entry>
        <entry key="bam">${bam}</entry>
        <entry key="bed">${bed}</entry>
</properties>
EOF

"""
}

process SORT_UNPAIRED_FASTQ {
tag "${key[0]} ${key[1]} N=${L.size()}"
afterScript "rm -f jeter.tsv.gz"
memory "3g"
input:
	val(meta)
	tuple val(key),val(L)
output:
	tuple val("${key[0]}"),val("${key[1]}"),path("${key[0]}.${key[1]}.tsv.gz"),emit:output
	path("version.xml"),emit:version
script:
"""
gunzip -c ${L.join(" ")} |\
	LC_ALL=C sort -T . -S ${task.memory.kilo} -t '\t' -k1,1 |\
	gzip --best > jeter.tsv.gz

mv jeter.tsv.gz "${key[0]}.${key[1]}.tsv.gz"


##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">sort unpaired linearized fastq on name</entry>
        <entry key="count">${L.size()}</entry>
        <entry key="sample">${key[0]}</entry>
        <entry key="side">${key[1]}</entry>
        <entry key="compressed.size">\$(ls -lah "${key[0]}.${key[1]}.tsv.gz")</entry>
</properties>
EOF
"""
}

process JOIN_UNPAIRED {
tag "${sample} ${R1} ${R2}"
afterScript "rm -f jeter.fq.gz jeter2.fq.gz"
memory "3g"
input:
	val(meta)
	tuple val(sample),val(R1),val(R2)
output:
	tuple val(sample),path("${sample}.R1R2.fq.gz"),emit:output
	tuple val(sample),path("${sample}.R0.fq.gz"),emit:single
	path("version.xml"),emit:version
script:
"""


LC_ALL=C join -t '\t' -1 1 -2 1 -o '1.1,1.2,1.3,1.4,2.1,2.2,2.3,2.4' \
	<(gunzip -c "${R1}")  \
	<(gunzip -c "${R2}") |\
	tr "\t" "\\n" |\
	gzip --best > jeter.fq.gz

LC_ALL=C join -t '\t' -1 1 -2 1 -v 1 -v2 \
	<(gunzip -c "${R1}")  \
	<(gunzip -c "${R2}") |\
	tr "\t" "\\n" |\
	gzip --best > jeter2.fq.gz


mv jeter.fq.gz  "${sample}.R1R2.fq.gz"
mv jeter2.fq.gz "${sample}.R0.fq.gz"

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">merge unpaired into interleaved and single</entry>
        <entry key="compressed.R1R2">\$(ls -lah "${sample}.R1R2.fq.gz")</entry>
        <entry key="compressed.R0">\$(ls -lah "${sample}.R0.fq.gz")</entry>
</properties>
EOF
"""
}
