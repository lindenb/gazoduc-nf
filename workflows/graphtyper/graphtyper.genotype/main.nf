/*

Copyright (c) 2026 Pierre Lindenbaum

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
include {JVARKIT_BAM_RENAME_CONTIGS               } from '../../../modules/jvarkit/bamrenamechr'


if(params.help) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}


Map assertKeyExists(final Map hash,final String key) {
    if(!hash.containsKey(key)) throw new IllegalArgumentException("no key ${key}'in ${hash}");
    return hash;
}

Map assertKeyExistsAndNotEmpty(final Map hash,final String key) {
    assertKeyExists(hash,key);
    def value = hash.get(key);
    if(value.isEmpty()) throw new IllegalArgumentException("empty ${key}'in ${hash}");
    return hash;
}

Map assertKeyMatchRegex(final Map hash,final String key,final String regex) {
    assertKeyExists(hash,key);
    def value = hash.get(key);
    if(!value.matches(regex)) throw new IllegalArgumentException(" ${key}'in ${hash} doesn't match regex '${regex}'.");
    return hash;
}

workflow {

	versions = Channel.empty()
	def genome = Channel.of( file(params.fasta), file(params.fai), file(params.dict) ).collect()	

  	bams = Channel.fromPath(params.samplesheet)
			.splitCsv(header:true,sep:',')
			.map{assertKeyMatchRegex(it,"sample","^[A-Za-z_0-9\\.\\-]+\$")}
			.map{assertKeyMatchRegex(it,"bam","^\\S+\\.(bam|cram)\$")}
			.map{
				if(it.containsKey("bai")) return it;
				if(it.bam.endsWith(".cram")) return it.plus(bai : it.bam+".crai");
				return it.plus(bai:it.bam+".bai");
			}
			.map{
				if(it.containsKey("fasta")) return it;
				return it.plus(fasta:params.fasta);
			}
			.map{
				if(it.containsKey("fai")) return it;
				return it.plus(fai:it.fasta+".fai");
			}
			.map{
				if(it.containsKey("dict")) return it;
				return it.plus(dict: it.fasta.replaceAll("\\.(fasta|fa|fna)\$",".dict"));
			}
			.map{
				if(it.containsKey("depth")) return it;
				return it.plus(depth:"-1");
			}
			.branch {
				ok_ref: it.fasta.equals(params.fasta)
				bad_ref: true
				}

	def bed_in= [];
	if(params.bed!=null) bed_in=file(params.bed)

	intervals_ch = MAKE_INTERVALS(genome,  bed_in)


	JVARKIT_BAM_RENAME_CONTIGS(
		[[id:"ref"],file(params.dict)],
		MAKE_INTERVALS.out.bed.map{[[id:"bed"],it]},
		bams.bad_ref.map{[[id:it.sample],file(it.bam),file(it.bai),file(it.fasta),file(it.fai),file(it.dict)]}
		)
	versions =versions.mix(JVARKIT_BAM_RENAME_CONTIGS.out.versions)

	bams_renamed = JVARKIT_BAM_RENAME_CONTIGS.out.bam
		.map{[it[0].id,it[1],it[2]]}
	



	bams2 = bams.ok_ref.branch{v->
			has_depth : v.containsKey("depth") && !v.depth.isEmpty() && !v.depth.equals(".") && (v.depth as int)>0
			no_depth : true
		}


	mosdepth_ch = MOSDEPTH(genome, intervals_ch.bed, bams2.no_depth.map{[it.sample,it.bam,it.bai]}.mix(bams_renamed) )

	ch1 = mosdepth_ch.output.splitText().
                map{[it[1],it[2],it[3],it[0].trim()]}.
		mix(bams2.has_depth.map{[it.sample,it.bam,it.bai,it.depth]})

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

	SAVE_FAILED(level6.failed.map{it.join("\t")}.collect())
	
        MERGE(merge_ch.map{T->{
			String c = T[0];
			if(!c.matches("(chr)?[0-9XY]+")) {
				c="others";
				}
			return [c,T[1]];
			}}.groupTuple())

        level6.failed.view{"${it} cannot be called"}
	
	}


runOnComplete(workflow);


process MAKE_INTERVALS {
label "process_single"
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
	def args1 = task.ext.args1 ?:" -w 50000 -s 49990"
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

cut -f1,2,3 TMP/jeter.bed |\\
	sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter2.bed
mv TMP/jeter2.bed TMP/jeter.bed


if ${bed?true:false}
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
	bedtools makewindows ${args1} -b - > TMP/windows.bed
test -s TMP/windows.bed
mv TMP/windows.bed ./
"""
}


process MOSDEPTH {
tag "${sample}"
label "process_single"
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
label "process_single"
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
set +o pipefail

samtools view -F 3844 -T "${fasta}" "${bam}" |\\
	cut -f 10 |\\
	head -n "${n_reads}" |\\
	awk -F '\t' 'BEGIN{T=0.0;N=0;} {N++;T+=length(\$1)} END{print (N==0?100:(T/N));}' > ${sample}.len.txt
"""
}

process DIGEST_DEPTH {
label "process_single"
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
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
        tuple val(contig),val(L)
output:
        path("${contig}.merged.bcf")
        path("${contig}.merged.bcf.csi")
script:
	def args = "--remove-duplicates"
"""
mkdir -p TMP
set -x
cat << EOF >  TMP/jeter.list
${L.join("\n")}
EOF

SQRT=`awk 'END{X=NR;z=sqrt(X); print (z==int(z)?z:int(z)+1);}' "TMP/jeter.list"`
split -a 9 --additional-suffix=.list --lines=\${SQRT} TMP/jeter.list TMP/chunck.


find TMP/ -type f -name "chunck*.list" | while read F
do
                bcftools concat --threads ${task.cpus} ${args} -a -O b --file-list "\${F}" -o "\${F}.bcf" 
		bcftools index --threads ${task.cpus} -f "\${F}.bcf"
                echo "\${F}.bcf" >> TMP/jeter2.list
done

bcftools concat ${args} --threads ${task.cpus} -a -O b9 --file-list TMP/jeter2.list -o "TMP/jeter.bcf" 
bcftools index --threads ${task.cpus} -f "TMP/jeter.bcf"

mv TMP/jeter.bcf ${contig}.merged.bcf
mv TMP/jeter.bcf.csi ${contig}.merged.bcf.csi
"""
}

process SAVE_FAILED {
executor "local"
input:
	val(L)
output:
	path("failed.txt"),emit:output
script:
"""
cat << EOF | sort > failed.txt
${L.join("\n")}
EOF
"""
}
