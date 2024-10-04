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

include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { HC_COMBINE1 } from '../modules/gatk/gatk4.combine.gvcfs.01.nf'
include { HC_COMBINE2 } from '../modules/gatk/gatk4.combine.gvcfs.02.nf'
include { HC_GENOTYPE } from '../modules/gatk/gatk4.genotype.gvcfs.01.nf'

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run my_pipeline --input input_file.csv")
   exit 0
}
// validate parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)


List makeSQRT(def L1) {
	def key = L1.get(0);
	def L = L1.get(1).sort();
	int n = (int)Math.ceil(Math.sqrt(L.size()));
	if(n<25) n=25;
	def returnList = [];
	def currList = [];
	int i=0;
	for(;;) {
		if(i<L.size()) currList.add(L.get(i));
		if(i==L.size() || currList.size()==n) {
			if(!currList.isEmpty()) returnList.add([key,currList]);
			if(i==L.size()) break;
			currList=[];
			}
		i++;
		}
	return returnList;
	}


workflow {

	beds_ch = Channel.fromPath(params.beds).
		splitText().
		map{file(it.trim())}




	genome_ch = [file(params.fasta) , file(params.fasta+".fai"), file(""+file(params.fasta).getParent()+"/"+file(params.fasta).getBaseName()+".dict" ) ]	
	if(params.dbsnp.equals("NO_FILE"))
		{
		dbsnp_ch = [file("NO_DBSNP"), file("NO_DBSNP_TBI")]
		}
	else
		{
		dbsnp_ch = [file(params.dbsnp), file(params.dbsnp+".tbi")]
		}

	hc_ch = Channel.empty()
        if(!params.gvcfs.equals("NO_FILE")) {
		gvcfs_ch = Channel.fromPath(params.gvcfs).
			splitText().
			map{it.trim()}.
			map{[it,it+".tbi"]}.
			map{[file(it[0],checkIfExists:true),file(it[1],checkIfExists:true)]}
		hc_ch = hc_ch.mix(  beds_ch.combine(gvcfs_ch) )
		}

	if(!params.bams.equals("NO_FILE")) {
		bams_ch = Channel.fromPath(params.bams).
			splitText().
			map{it.trim()}.
			map{[it, it.endsWith(".cram")?it+".crai":it+".bai"]}.
			map{[file(it[0],checkIfExists:true),file(it[1],checkIfExists:true)]}

		hc_ch = hc_ch.mix( HC_BAM_BED( genome_ch, dbsnp_ch, beds_ch.combine(bams_ch)) )
		}
	

	combine1_input_ch = hc_ch.output.map{[it[0].toRealPath(),[it[1],it[2]]]}.
		groupTuple().
		flatMap{makeSQRT(it)}.
		map{[it[0],it[1].flatten()]}

	hc1 = HC_COMBINE1(genome_ch, dbsnp_ch, combine1_input_ch )

	combine2_input_ch = hc1.output.
		map{ [it[0].toRealPath(), it[1],it[2] ]}. /* duplicate bed file so we can extract distinct contig in first , keep bed content in second */
		groupTuple(). /* group by contig/bed */
		map{[it[0], it[1].flatten().plus(it[2].flatten())]}

	splitbed_ch = SPLIT_BED(combine2_input_ch)

	combine2_input_ch = splitbed_ch.output.
		flatMap{T->T[0].collect{X->[X,T[1]]} }

	hc2 = HC_COMBINE2(genome_ch, dbsnp_ch, combine2_input_ch )
	hg = HC_GENOTYPE(genome_ch, dbsnp_ch, hc2.output )
	concat_ch = VCF_CONCAT(hg.output.flatten().collect())
	stats_ch = BCFTOOLS_STATS(genome_ch,concat_ch.output)
	mqc_ch = MULTIQC(file(params.sample2population), stats_ch.output)
	}

process HC_BAM_BED {
tag "${bed.name} ${bam.name}"
label "process_quick"
afterScript "rm -rf TMP"
errorStrategy "retry"
maxRetries 2
cpus 2
input:
	tuple path(fasta),path(fai),path(dict)
	tuple path(dbsnp),path(dbsnp_tbi)
	tuple path(bed),path(bam),path(bai)
output:
	tuple path(bed),path("${bam.getBaseName()}.g.vcf.gz"),path("${bam.getBaseName()}.g.vcf.gz.tbi"),emit:output
script:
"""
hostname 1>&2
module load gatk/0.0.0
mkdir -p TMP

   gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" HaplotypeCaller \\
     -L "${bed}" \\
     -R "${fasta}" \\
     -I "${bam}" \\
     -ERC GVCF \\
     ${dbsnp.name.equals("NO_DBSNP")?"":"--dbsnp ${dbsnp}"} \\
     ${(params.mapq as Integer)<1?"":" --minimum-mapping-quality "+params.mapq} \\
     -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \\
     -O "TMP/jeter.g.vcf.gz"

mv TMP/jeter.g.vcf.gz "${bam.getBaseName()}.g.vcf.gz"
mv TMP/jeter.g.vcf.gz.tbi "${bam.getBaseName()}.g.vcf.gz.tbi"
"""
}


process SPLIT_BED {
tag "${bed.name}"
executor "local"
maxForks 10
input:
	tuple path(bed),path(vcfs)
output:
	tuple path("BEDS/*.bed"),path(vcfs),emit:output
script:
	def f = bed.getBaseName();
	def size="1E6";
"""
mkdir -p BEDS
awk -F '\t' 'BEGIN{N=1;T=0;f=sprintf("BEDS/${f}.%d.bed",N);} {print \$0 >> f; T+=int(\$3)-int(\$2); if(T>${size}) {close(f);T=0;N++;f=sprintf("BEDS/${f}.%d.bed",N);}}' '${bed}'
"""
}


process VCF_CONCAT {
label "process_high"
cpus 20
time "24h"
afterScript "rm -rf TMP"
input:
	path("VCFS/*")
output:
	tuple path("output.bcf"),path("output.bcf.csi"),emit:output
script:
	
"""
hostname 1>&2
module load bcftools
set -o pipefail
mkdir -p TMP
set -x
	
find ./VCFS -name "*.vcf.gz" > TMP/jeter.list

bcftools concat --threads ${task.cpus} --allow-overlaps --rm-dups exact --file-list TMP/jeter.list -O b -o  TMP/jeter.bcf
bcftools sort -T TMP/sort. --max-mem "${(task.memory.giga / 2 ) as int }G" -O b9 -o TMP/output.bcf TMP/jeter.bcf
bcftools index --threads ${task.cpus} -f TMP/output.bcf
mv -v TMP/output.bcf ./
mv -v TMP/output.bcf.csi ./
"""
}

process BCFTOOLS_STATS {
label "process_medium"
cpus 10
input:
        tuple path(fasta),path(fai),path(dict)
	tuple path(vcf),path(vcfidx)
output:
	path("stats.txt"),emit:output	
script:
	
"""
hostname 1>&2
module load  bcftools
mkdir -p TMP

bcftools stats --threads ${task.cpus} --samples - --fasta-ref "${fasta}" "${vcf}" > "stats.txt"

"""
}

process MULTIQC {
label "process_medium"
input:
	path(sample2pop)
	path(stats)
output:
	path("multiqc.zip"),emit:output
script:
	def prefix = params.prefix
"""
hostname 1>&2
module load jvarkit multiqc
mkdir -p TMP
echo "${stats}" > TMP/jeter.list

export LC_ALL=en_US.utf8

if ${!sample2pop.name.equals("NO_FILE")}
then

mkdir -p TMP/OUT2 TMP/multiqc

multiqc --outdir "TMP/multiqc" --force --file-list TMP/jeter.list --no-ansi

java -jar \${JVARKIT_DIST}/jvarkit.jar multiqcpostproc --sample2collection "${sample2pop}" -o TMP/OUT2 TMP/multiqc/multiqc_data

find TMP/OUT2 -type f -name "*.json" >> TMP/jeter.list

fi

		
	mkdir -p "${prefix}multiqc"
	multiqc  --filename  "${prefix}multiqc_report.html" --no-ansi \
			--title "${prefix}GATK4 Calling"  \
			--comment "calling with gatk4 ${prefix}.gatk"  \
			--force \
			--outdir "${prefix}multiqc" \
			--file-list TMP/jeter.list
		
	rm -f multiqc.zip
	zip -9 -r "multiqc.zip" "${prefix}multiqc"

"""
}
