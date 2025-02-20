/*

Copyright (c) 2025 Pierre Lindenbaum

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

//include {GRAPHTYPER_GENOTYPE_BAMS_01} from '../../../subworkflows/graphtyper/graphtyper.genotype.bams.01.nf'
include {dumpParams;runOnComplete} from '../../../modules/utils/functions.nf'
include {DIVIDE_AND_CONQUER as DIVIDE_AND_CONQUER1} from "./sub.nf"
include {DIVIDE_AND_CONQUER as DIVIDE_AND_CONQUER2} from "./sub.nf"
include {DIVIDE_AND_CONQUER as DIVIDE_AND_CONQUER3} from "./sub.nf"
include {DIVIDE_AND_CONQUER as DIVIDE_AND_CONQUER4} from "./sub.nf"
include {DIVIDE_AND_CONQUER as DIVIDE_AND_CONQUER5} from "./sub.nf"
include {DIVIDE_AND_CONQUER as DIVIDE_AND_CONQUER6} from "./sub.nf"


if(params.help) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}


workflow {


	def genome = Channel.of( file(params.fasta), file(params.fai), file(params.dict) ).collect()	

	bams_ch = Channel.fromPath(params.samplesheet).
		splitCsv(header:true,sep:'\t').
		branch{v->
			has_depth : v.containsKey("depth") && !v.depth.isEmpty() && !v.depth.equals(".")
			no_depth : true
		}


	intervals_ch = MAKE_INTERVALS(genome, file(params.bed))
	mosdepth_ch = MOSDEPTH(genome, intervals_ch.bed, bams_ch.no_depth.map{[it.sample,it.bam,it.bai]} )

	ch1 = mosdepth_ch.output.splitText().
                map{[it[1],it[2],it[3],it[0].trim()]}.
		mix(bams_ch.has_depth.map{[it.sample,it.bam,it.bai,it.depth]})

	readlen_ch = READ_LENGTH( genome, ch1.map{[it[0],it[1]]} )


	ch2 = readlen_ch.splitText().map{[it[1],it[0].trim()]}.join(ch1)

	samplesheet_ch = DIGEST_DEPTH(ch2.map{it.join("\t")}.collect())


	
        intervals_ch  = intervals_ch.output.
                splitCsv(header:false,sep:'\t').
                map{[it[0],""+((it[1] as int)+1),it[2]]}


        merge_ch = Channel.empty()

        level1 = DIVIDE_AND_CONQUER1(1, samplesheet_ch.output , intervals_ch )
        merge_ch = merge_ch.mix(level1.ok)

        level2 = DIVIDE_AND_CONQUER2(2, samplesheet_ch.output , level1.failed)
        merge_ch = merge_ch.mix(level2.ok)

        level3 = DIVIDE_AND_CONQUER3(3, samplesheet_ch.output , level2.failed)
        merge_ch = merge_ch.mix(level3.ok)

        level4 = DIVIDE_AND_CONQUER4(4, samplesheet_ch.output , level3.failed)
        merge_ch = merge_ch.mix(level4.ok)

        level5 = DIVIDE_AND_CONQUER5(5, samplesheet_ch.output , level4.failed)
        merge_ch = merge_ch.mix(level5.ok)

        level6 = DIVIDE_AND_CONQUER6(6, samplesheet_ch.output , level5.failed)
        merge_ch = merge_ch.mix(level6.ok)


        MERGE(merge_ch.groupTuple())

        level6.failed.view{"${it} cannot be called"}

	}


runOnComplete(workflow);


process MAKE_INTERVALS {
label "process_short"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	path(genome)
	path(bed)
output:
	path("windows.bed"),emit:output
	path("all.bed"),emit:bed
script:
	def fasta = genome.find{it.name.endsWith("a")}
	def regex= '^(chr)?[0-9XY]+\$'
"""
mkdir -p TMP
export LC_ALL=C

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" ScatterIntervalsByNs \\
	    --REFERENCE "${fasta}" \\
	    --MAX_TO_MERGE "1000" \\
	    --OUTPUT "TMP/jeter.interval_list" \\
	    --OUTPUT_TYPE "ACGT"

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" IntervalListToBed \\
	    --INPUT "TMP/jeter.interval_list" \\
	    --OUTPUT "TMP/jeter.bed" \\
	    --SORT true

if ${bed.name.contains(".")}
then
	cut -f1,2,3 "${bed}" | sort -T TMP -t '\t' -k1,1 -k2,2n | bedtools merge > TMP/jeter2.bed
	bedtools intersect -a TMP/jeter.bed -b TMP/jeter2.bed > TMP/jeter3.bed
	mv TMP/jeter3.bed TMP/jeter.bed
fi

awk -F '\t' '(\$1 ~ /${regex}/ )' TMP/jeter.bed |\\
	cut -f1,2,3 |\\
	sort -T TMP -t '\t' -k1,1 -k2,2n > all.bed
	
cut -f1,2,3 all.bed |\\
	sort -T TMP -t '\t' -k1,1 -k2,2n  |\\
	bedtools makewindows -w 50000 -s 49990 -b - > TMP/windows.bed
test -s TMP/windows.bed
mv TMP/windows.bed ./
"""
}


process MOSDEPTH {
tag "${sample}"
label "process_quick"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/mosdepth.yml"
input:
	path(genome)
	path(bed)
	tuple val(sample),path(bam),path(bai)
output:
	tuple path("${sample}.cov.txt"),val(sample),path(bam),path(bai), emit:output
script:
	def fasta = genome.find{it.name.endsWith("a")}
	def mapq = params.mapq
"""
mkdir -p TMP

# bed for autosomes
awk -F '\t' '(\$1 ~/^(chr)?[0-9]+\$/)' '${bed}' > TMP/jeter.bed
test -s TMP/jeter.bed

mosdepth  \\
 	-t ${task.cpus} \\
	--by TMP/jeter.bed \\
	--no-per-base \\
	--fasta "${fasta}" \\
	--mapq ${mapq} \\
	TMP/output \\
	${bam}

awk -F '\t' '(\$1=="total_region") {print \$4}' TMP/output.mosdepth.summary.txt > "${sample}.cov.txt"
test "${sample}.cov.txt"
"""
}



process READ_LENGTH {
tag "${sample}"
label "process_quick"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	path(genome)
	tuple val(sample),path(bam)
output:
	tuple path("${sample}.len.txt"),val(sample)
script:
	def fasta = genome.find{it.name.endsWith("a")}
	def n_reads = 1000
"""

mkdir -p TMP

samtools view -F 3844 -T "${fasta}" "${bam}" |\\
	cut -f 10 |\\
	head -n "${n_reads}" |\\
	awk -F '\t' 'BEGIN{T=0.0;N=0;} {N++;T+=length(\$1)} END{print (N==0?100:(T/N));}' > ${sample}.len.txt
"""
}

process DIGEST_DEPTH {
label "process_quick"
input:
	val(L)
output:
	path("samplesheet.tsv"),emit:output
script:
"""
cat << EOF > jeter.tsv
${L.join("\n")}
EOF

sort -T . -t '\t' -k1,1 jeter.tsv > samplesheet.tsv
test -s samplesheet.tsv
cut -f1 samplesheet.tsv | sort | uniq -d > dups.txt
test ! -s dups.txt
"""
}


process RDF {
input:
	val(L)
output:
	path("digest.ttl"),emit:output
script:
"""

cat << EOF | awk -F '\t' '{printf("samples:%s rdf:type foaf:Person;\\n\tu:Library libraries:\\"%s\".\\n",\$1,\$1); }' > digest.ttl
${L.join("\n")}
EOF

"""
}


process MERGE {
tag "${contig} N=${L.size()}"
label "process_medium"
afterScript "rm -rf TMP"
input:
        tuple val(contig),val(L)
output:
        path("${contig}.merged.bcf")
        path("${contig}.merged.bcf.csi")
script:
"""
module load bcftools
mkdir -p TMP
set -x
cat << EOF >  TMP/jeter.list
${L.join("\n")}
EOF

SQRT=`awk 'END{X=NR;z=sqrt(X); print (z==int(z)?z:int(z)+1);}' "TMP/jeter.list"`
split -a 9 --additional-suffix=.list --lines=\${SQRT} TMP/jeter.list TMP/chunck.


find TMP/ -type f -name "chunck*.list" | while read F
do
                bcftools concat --write-index --threads ${task.cpus} -a -O b --file-list "\${F}" -o "\${F}.bcf" 
                echo "\${F}.bcf" >> TMP/jeter2.list
done

bcftools concat --write-index --threads ${task.cpus} -a -O b9 --file-list TMP/jeter2.list -o "TMP/jeter.bcf" 

mv TMP/jeter.bcf ${contig}.merged.bcf
mv TMP/jeter.bcf.csi ${contig}.merged.bcf.csi
"""
}
