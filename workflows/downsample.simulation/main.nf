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


include {assertKeyExistsAndNotEmpty          } from '../../modules/utils/functions.nf'
include {PREPARE_REFERENCE                   } from '../../subworkflows/samtools/prepare.ref'
include {META_TO_PED                         } from '../../subworkflows/pedigree/meta2ped'
include {MULTIQC                             } from '../../subworkflows/multiqc'
include {META_TO_BAMS                        } from '../../subworkflows/samtools/meta2bams1'
include {runOnComplete                       } from '../../modules/utils/functions.nf'
include {MOSDEPTH as MOSDEPTH1               } from '../../modules/mosdepth'
include {MOSDEPTH as MOSDEPTH2               } from '../../modules/mosdepth'
include {BEDTOOLS_INTERSECT                  } from '../../modules/bedtools/intersect'
include {DOWNSAMPLE                          } from '../../modules/samtools/downsample'
include {SAMTOOLS_VIEW                       } from '../../modules/samtools/view'
include {HAPLOTYPECALLER                     } from '../../modules/gatk/hapcaller1'
include {GENOTYPEGVCFS                       } from '../../modules/gatk/genotypegvcfs'
include {COMBINEGVCFS                        } from '../../modules/gatk/combinegvcfs'
include {GHOSTSCRIPT_MERGE                   } from '../../modules/gs/merge'
include {PREPARE_USER_BED                    } from '../../subworkflows/bedtools/prepare.user.bed'

workflow {
		if(params.samplesheet==null) {
			throw new IllegalArgumentException("--samplesheet not defined");
			}
		if(params.fasta==null) {
			throw new IllegalArgumentException("--fasta not defined");
			}

        versions = Channel.empty()
        multiqc  = Channel.empty()
		each_depths = Channel.of(params.depths.split("[,]"))
			.flatMap()
			.map{T->(T as int)}
			.filter{it>0}
		
        def hash_ref= [
                id: file(params.fasta).baseName,
                name: file(params.fasta).baseName,
                ucsc_name: (params.ucsc_name?:"undefined")
                ]
        PREPARE_REFERENCE(
			hash_ref.plus(skip_complement:true),
			Channel.of( [ hash_ref, file(params.fasta)])
			)
        versions = versions.mix(PREPARE_REFERENCE.out.versions)
		fasta = PREPARE_REFERENCE.out.fasta.first()
		fai = PREPARE_REFERENCE.out.fai.first()
		dict = PREPARE_REFERENCE.out.dict.first()

        META_TO_BAMS(
			hash_ref,
			fasta,
			fai,
			Channel.fromPath(params.samplesheet)
				.splitCsv(header:true, sep:',')
			)
        versions = versions.mix(META_TO_BAMS.out.versions)
	    bams_ch = META_TO_BAMS.out.bams

		if(params.bed==null) {
			bed = Channel.empty()
			}
		else
			{
			bed = Channel.of([[id:file(params.bed).baseName],file(params.bed)])
			}

		PREPARE_USER_BED(
			hash_ref,
			fasta,
			fai,
			dict,
			PREPARE_REFERENCE.out.scatter_bed.first(),
			bed
			)
		versions = versions.mix(PREPARE_USER_BED.out.versions)
		multiqc = multiqc.mix(PREPARE_USER_BED.out.multiqc)
		bed = PREPARE_USER_BED.out.bed.first()
  		
  		/*  will force the original BAM to be filtered as the downsampled one */
	    		
		if(params.use_original_bam!=true) {

			SAMTOOLS_VIEW(
				fasta,
				fai,
				bed,
				[[id:"no_reads"],[]],
				bams_ch
				);
			versions = versions.mix(SAMTOOLS_VIEW.out.versions)
			
			bams_ch = SAMTOOLS_VIEW.out.bam
			}
  		
  		MOSDEPTH1(
			fasta,
			fai,
			bams_ch
				.combine(bed)
				.map{meta,bam,bai,meta2,bed->[meta,bam,bai,bed]}
			)
  		versions = versions.mix(MOSDEPTH1.out.versions)
		multiqc = multiqc
				.mix(MOSDEPTH1.out.summary_txt)
				.mix(MOSDEPTH1.out.global_txt)
				.mix(MOSDEPTH1.out.global_txt)
  
  		ch1 = bams_ch.map{meta,bam,bai->[meta.id,meta,bam,bai]}
  		ch2 = MOSDEPTH1.out.summary_txt.map{meta,summary->[meta.id,summary]}
  		
  		
  		ch1 = ch1.join(ch2)
  			.map{meta_id,meta,bam,bai,summary->[meta,bam,bai,summary]}
  			.combine(each_depths)
  			.map{meta,bam,bai,summary,dp->[meta.plus(depth:dp),bam,bai,summary]}
  			
  		DOWNSAMPLE(
  			fasta,
  			fai,
  			bed,
  			ch1
  			)
  	   versions = versions.mix(DOWNSAMPLE.out.versions)
  		
       MOSDEPTH2(
		fasta,
		fai,
		DOWNSAMPLE.out.bam
			.combine(bed)
			.map{meta,bam,bai,meta2,bed->[meta,bam,bai,bed]}
		)
  	   versions = versions.mix(MOSDEPTH2.out.versions)
	   multiqc = multiqc
				.mix(MOSDEPTH2.out.summary_txt)
				.mix(MOSDEPTH2.out.global_txt)
				.mix(MOSDEPTH2.out.global_txt)
	
		bams_ch= bams_ch
		   		.mix(DOWNSAMPLE.out.bam)
		   		.map{meta,bam,bai->{[
		   			meta.plus(collection: (meta.depth?"DP"+meta.depth:"RAW"), original_sample: meta.id), // ADD COLLECTION, SRC SAMPLE
		   			bam,
		   			bai
		   			]}}
		   		.map{meta,bam,bai->{[
		   			meta.plus(id: meta.id+(meta.depth?".DP"+meta.depth:"")), // MODIFY ID for GATK
		   			bam,
		   			bai
		   			]}}
		 
		META_TO_PED(hash_ref, bams_ch.map{it[0]})
        versions = versions.mix(META_TO_PED.out.versions)
		 
	
	   HAPLOTYPECALLER(
		   fasta,
		   fai,
		   dict,
		   [[id:"noref"],[]],
		   bams_ch
		   		.combine(bed)
		   		.map{meta,bam,bai,meta2,bed->[meta,bam,bai,bed]}
		   )
	  versions = versions.mix(HAPLOTYPECALLER.out.versions)
	  
	  

	  gvcf_ch = HAPLOTYPECALLER.out.gvcf
		    	.map{meta,vcf,tbi,bed->[meta.original_sample,[vcf,tbi],bed]}
		    	.groupTuple()
		    	.map{sn,vcf_files,beds->[[id:sn],vcf_files.flatten().sort(),beds[0]]}
		    	

	  COMBINEGVCFS(
			fasta,
			fai,
		    dict,
		    gvcf_ch
		    )
	  versions = versions.mix(COMBINEGVCFS.out.versions)
	  GENOTYPEGVCFS(
	  		fasta,
			fai,
		    dict,
		    [[id:"nodbsnp"],[],[]],
		    COMBINEGVCFS.out.gvcf
	  		)
	  versions = versions.mix(GENOTYPEGVCFS.out.versions)
	
	
	GENOTYPE_CONCORDANCE(GENOTYPEGVCFS.out.vcf)
	versions = versions.mix(GENOTYPE_CONCORDANCE.out.versions)

	

	PLOT(GENOTYPE_CONCORDANCE.out.output.map{it[1]}.collect().sort().map{[[id:"lowpass"],it]})
	versions = versions.mix(PLOT.out.versions)


	GHOSTSCRIPT_MERGE(PLOT.out.pdf.mix(GENOTYPE_CONCORDANCE.out.pdf).map{it[1]}.collect().sort().map{[[id:"lowpass"],it]})
	versions = versions.mix(GHOSTSCRIPT_MERGE.out.versions)
	
	MULTIQC(
            hash_ref.plus("id":"downsample"),
            META_TO_PED.out.sample2collection,
            versions,
			[[id:"no_mqc_config"],[]],
            multiqc
            )
	}

runOnComplete(workflow)



process GENOTYPE_CONCORDANCE {
   tag "${meta.id} N=${vcf.name}"
   conda "${moduleDir}/../../conda/bioinfo.02.yml"
   afterScript "rm -rf TMP"
   input:
        tuple val(meta),path(vcf),path(tbi),path(bed)
   output:
		tuple val(meta),path("*.concordances.txt"),emit:output
		tuple val(meta),path("*.GQ.pdf"),emit:pdf
        path("versions.yml"),emit:versions
   script:
	def jvm= task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"
   """
   hostname 1>&2
   set -o pipefail
   mkdir -p TMP

   bcftools query -l "${vcf}" | grep -v -x '${meta.id}' | while read S
   do

   gatk --java-options "${jvm}" GenotypeConcordance \\
	--TRUTH_VCF "${vcf}" \\
	--TRUTH_SAMPLE "${meta.id}" \\
	--CALL_VCF "${vcf}" \\
	--CALL_SAMPLE "\${S}" \\
	--O "\${S}"

    grep -E '^(SNP|INDEL)' "\${S}.genotype_concordance_summary_metrics"  |\\
	cut -f 1,13 |\\
	awk -vS=\$S '{printf("%s\t%s\\n", gsub(/.*\\.DP[0]*/,"\\\\1",S) ,\$0);}' >> TMP/concordances.txt
   done


   bcftools query -l "${vcf}" | while read S
   do
      bcftools view -O u --samples "\${S}" "${vcf}" |\\
		bcftools query -f '[%GQ\\n]' | sed 's/^\\.\$/-1/' | LC_ALL=C sort -n >> "TMP/\${S}.dist"
   done


echo -n "filenames <-c(" > TMP/jeter.R

find TMP -type f -name "*.dist" | sort | awk '{printf("\\"%s\\"\\n",\$0);}' |paste -s -d, >> TMP/jeter.R

cat << '__EOF__' >> TMP/jeter.R
)
n <- 1
colors <- rainbow(length(filenames),alpha=0.8)
pdf("${meta.id}.GQ.pdf")
list1 <- list()
list2 <- c()
maxy <- 1E-6

for(filename in filenames) {
        X <- data.frame(table(scan(filename, comment.char = ".")))
	maxy <- max(maxy, max(X\$Freq))
}

for(filename in filenames) {
	T1 <- table(scan(filename, comment.char = "."))
	X <- data.frame(T1)
	list2 <- append(list2,gsub("${meta.id}.DP","",gsub(".dist","",basename(filename))))

	if(n==1) {
		plot(x=as.integer(as.character(X\$Var1)),y=X\$Freq,type='l',main="${meta.id} : density of Genotype Quality",xlab="Genotype Quality",ylab="Count",col=colors[n],ylim=c(0,maxy),xlim=c(0,100))
	} else	{		
		lines(x=as.integer(as.character(X\$Var1)),y=X\$Freq,col=colors[n],ylim=c(0,maxy),xlim=c(0,100))
		}
	n <- n+1
}

legend("top", legend=list2,text.col=colors)
dev.off()
__EOF__

R --vanilla < TMP/jeter.R
mv TMP/concordances.txt "${meta.id}.concordances.txt"
###############################################################################
touch versions.yml
"""
stub:
"""
touch versions.yml ${meta.id}.concordances.txt ${meta.id}.GQ.pdf
"""
}

process PLOT {
   label "process_single"
   conda "${moduleDir}/../../conda/bioinfo.01.yml"
   afterScript "rm -rf TMP"
input:
	tuple val(meta),path("FILES/*")
output:
	tuple val(meta),path("concordance.pdf"),emit:pdf
	path("versions.yml"),emit:versions
script:
"""
hostname 1>&2
mkdir -p TMP

find FILES/ -name "*.txt" -exec cat '{}' ';' > TMP/jeter.tsv
test -s TMP/jeter.tsv

cat << '__EOF__' > TMP/jeter.R

T1 <- read.table("TMP/jeter.tsv",sep="\t",header=FALSE)
depths <- sort(as.integer(unique(T1\$V1)), decreasing = TRUE)
types <- unique(T1\$V2)

list1 <- list()
list2 <- c()
list3 <- c()

for(type in types) {
for(depth in depths) {

	T2 <-T1[as.integer(T1\$V1)==depth & T1\$V2==type,]\$V3
	list1[[length(list1)+1]] <- as.double(T2)
	list2 <- append(list2,depth)
	list3 <- append(list3,ifelse(type=="SNP","blue","yellow"))
}
}


pdf("TMP/concordance.pdf")

boxplot(list1,names=list2,main="concordance",xlab="depth",ylab="type",col=list3,las=2)

legend("topright",legend=c("SNP","INDEL"),fill=c("blue","yellow"))

dev.off()

__EOF__

R --vanilla < TMP/jeter.R

mv TMP/concordance.pdf ./
touch versions.yml
"""

stub:
"""
touch versions.yml "concordance.pdf"
"""
}

