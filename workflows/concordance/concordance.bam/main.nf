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
include { runOnComplete; dumpParams   } from '../../../modules/utils/functions.nf'
include { PREPARE_ONE_REFERENCE       } from '../../../subworkflows/samtools/prepare.one.ref'
include { GHOSTSCRIPT_MERGE           } from '../../../modules/gs/merge'

workflow {
	if(params.samplesheet==null) {
		log.warn("undefined --samplesheet");
		exit -1
		}
	if(params.fasta==null) {
		log.warn("undefined --fasta");
		exit -1
		}
	if(params.bed==null) {
		log.warn("undefined --bed");
		exit -1
		}
    if(params.gnomad==null) {
		log.warn("undefined --gnomad");
		exit -1
		}
   versions = Channel.empty()
   multiqc = Channel.empty()
   metadata = [
		id:"concordance"
		]

   PREPARE_ONE_REFERENCE(
   		[id:"x",skip_scatter:true],
   		Channel.of(params.fasta).map{f->file(f)}.map{f->[[id:f.baseName],file(f)]}.first()
   		)

   	bed_ch =	Channel.of(params.bed).map{f->file(f)}.map{f->[[id:f.baseName],file(f)]}.first()
    gnomad_ch =	Channel.of(params.gnomad).map{f->file(f)}.map{f->[[id:f.baseName],file(f),file(f+".tbi")]}.first()
   
   	bams = Channel.fromPath(params.samplesheet)
   		.splitCsv(header:true,sep:',')
   		.map{row->[
				[
					id:row.id,
					type:row.type,
					sample:row.sample
				],
				file(row.bam),
				file(row.bam+".crai"),
				]}
    GENOTYPE_TRUTH(
   		PREPARE_ONE_REFERENCE.out.fasta,
   		PREPARE_ONE_REFERENCE.out.fai,
   		PREPARE_ONE_REFERENCE.out.dict,
   		bed_ch,
        gnomad_ch,
   		bams.filter{meta,_bam,_bai->meta.type=="truth"}
   		)
    

   	CONCAT_SITES(
   		PREPARE_ONE_REFERENCE.out.fasta,
   		PREPARE_ONE_REFERENCE.out.fai,
   		GENOTYPE_TRUTH.out.vcf
   			.map{meta,vcf,tbi->vcf}
   			.collect()
   			.map{[[id:"merge"],it.sort()]},
   		)
   		

   	GENOTYPE_CALL(
   		PREPARE_ONE_REFERENCE.out.fasta,
   		PREPARE_ONE_REFERENCE.out.fai,
   		PREPARE_ONE_REFERENCE.out.dict,
   		CONCAT_SITES.out.vcf,
   		bams.filter{meta,_bam,_bai->meta.type!="truth"}
   		)

    BCFTOOLS_MERGE(
   		GENOTYPE_CALL.out.vcf
			.map{meta,vcf,tbi->[vcf,tbi]}
			.flatMap()
			.collect()
			.map{[[id:"merge"],it.sort()]}
   		)
    CONCORDANCE(
   		GENOTYPE_TRUTH.out.vcf,
   		BCFTOOLS_MERGE.out.vcf
   		)
	versions = versions.mix(CONCORDANCE.out.versions)

	PLOT_CONCORDANCE(CONCORDANCE.out.concordance)
	versions = versions.mix(PLOT_CONCORDANCE.out.versions)


   	GHOSTSCRIPT_MERGE(
		PLOT_CONCORDANCE.out.pdf
			.map{meta,f->f}
			//.map{[it.name,it]} //PREVENT FILE NAME COLLISTION
			//.groupTuple()
			//.map{f,L->L[0]}
			.collect()
			.map{files->[ [id:"plot"],files.sort()]}
		)
   	 versions = versions.mix(GHOSTSCRIPT_MERGE.out.versions)

   	 MERGE_CONCORDANCES( 
		CONCORDANCE.out.concordance
			.map{meta,f->f}
			.collect()
			.map{files->[[id:"concordance"],files.sort()]}
		)
	versions = versions.mix(MERGE_CONCORDANCES.out.versions)

	README(
		bed_ch ,
		gnomad_ch ,
		GENOTYPE_TRUTH.out.vcf.first(),
		CONCAT_SITES.out.vcf ,
		GENOTYPE_CALL.out.vcf.first() ,
		BCFTOOLS_MERGE.out.vcf ,
		bams.filter{meta,_bam,_bai->meta.type=="truth"}
			.map{meta,bam,bai->bam}
			.take(10)
			.collect()
			.map{f->[[id:"truth"],f.sort()]} ,
		bams.filter{meta,_bam,_bai->meta.type!="truth"}
			.map{meta,bam,bai->bam}
			.take(10)
			.collect()
			.map{f->[[id:"cases"],f.sort()]} ,
		CONCORDANCE.out.concordance.first()
	)

}

process README {
executor "local"
input:
	tuple val(meta1),path(bed)
	tuple val(meta2),path(gnomad),path(gnomad_tbi)
	tuple val(meta2a),path(truth_vcf),path(tbi1)
	tuple val(meta2b),path(concat_vcf),path(concat_vcf_tbi)
	tuple val(meta2c),path(call_vcf),path(tbi2)
	tuple val(meta2d),path(merge_vcf),path(tbi3)
	tuple val(meta3),val(truth)
	tuple val(meta4),val(calls)
	tuple val(meta5),path(concordance)
output:
	path("README.md")
script:
	def antislash = '`'
	def antislash3 = '```'
"""
set +o pipefail

cat << '__EOF__' > README.md
# BAM concordance
__EOF__

############################################################

if ${bed?true:false}
then

cat << _'_EOF__' >> README.md

## Input BED

Input bed was ${antislash}${bed.toRealPath()}${antislash}.

${antislash3}
__EOF__

head "${bed}" >> README.md

cat << '__EOF__' >> README.md
(...)
${antislash3}
__EOF__

fi


############################################################


cat << '__EOF__' >> README.md

## Input GNOMAD

Gnomad was ${antislash}${gnomad.toRealPath()}${antislash}.

__EOF__


############################################################
cat << '__EOF__' >> README.md

## TRUTH Samples:

example:
${antislash3}
__EOF__

cat << EOF >> README.md
${truth.join("\n")}
EOF

cat << '__EOF__' >> README.md
(...)
${antislash3}
__EOF__

############################################################
cat << '__EOF__' >> README.md

## CALL Samples:

example:
${antislash3}
__EOF__

cat << 'EOF' >> README.md
${calls.join("\n")}
EOF

cat << '__EOF__' >> README.md
(...)
${antislash3}
__EOF__
############################################################


cat << '__EOF__' >> README.md

## Workflow

'Truth' samples are genotyped using GATK Haplotype caller in ${bed.name}.
We keep the snp diallelic variants with het genotypes  with a depth > 4 and depth <= 30 and MQ>=55.
Variants are discarded if they're in gnomad with AF/NFE > 0.01


For example:
${antislash3}
__EOF__


gunzip -c "${truth_vcf}" | grep -v "^##" | head >> README.md


cat << '__EOF__' >> README.md
(...)
${antislash3}

All 'Truth' vcfs are combined to create a site-only VCF site 'truth.vcf' containing
all the variants that have been found in all 'Truth' VCF.


${antislash3}
__EOF__

gunzip -c "${concat_vcf}" | grep -v "^##" | head >> README.md


cat << '__EOF__' >> README.md
(...)
${antislash3}


The 'truth.vcf' is used to **genotype** the **CALL** bams.
For example:
${antislash3}
__EOF__

gunzip -c "${call_vcf}" | grep -v "^##" | head >> README.md


cat << '__EOF__' >> README.md
(...)
${antislash3}

All 'Call' vcf are merged into one 'call.vcf' that contains all the genotypes.

${antislash3}
__EOF__

gunzip -c "${merge_vcf}" | grep -v "^##" | head | cut -f 1-20 | sed 's/\$/ (...)/' >> README.md

cat << '__EOF__' >> README.md
(...)
${antislash3}

Each 'Truth' vcf is challenged with 'bcftools gtcheck' against the 'call.vcf'
The output looks like:
${antislash3}
__EOF__

head ${concordance} >> README.md

cat << '__EOF__' >> README.md
${antislash3}

For each 'TRUTH' samples we plot the number of informative sites vs the concordance of the CALL samples.
__EOF__

"""
}


process GENOTYPE_TRUTH {
tag "${meta.sample}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta4),path(bed)
    tuple val(meta5),path(gnomad),path(gnomad_tbi)
	tuple val(meta ),path(bam),path(bai)
output:
	tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),optional:true,emit:vcf
script:
   def jvm = task.ext.jvm?:"-XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -DGATK_STACKTRACE_ON_USER_EXCEPTION=true"
   def max_af = (task.ext.max_af?:0.01) as double
   def prefix = task.ext.prefix ?:"${meta.sample}.${meta.type}."
"""
set -x
mkdir -p TMP/isec1 TMP/isec2

gatk --java-options "${jvm}" HaplotypeCaller \\
    -L "${bed}" \\
    -R "${fasta}" \\
    -I "${bam}" \\
    -O "TMP/jeter.vcf.gz" \\

# select good one
bcftools view -m 2 -M 2 -i 'AC=1 && INFO/DP>=4 && INFO/DP<=30 && INFO/MQ>=55' --type snps -O z -o TMP/jeter2.vcf.gz TMP/jeter.vcf.gz
bcftools index -f -t TMP/jeter2.vcf.gz
mv TMP/jeter2.vcf.gz  TMP/jeter.vcf.gz
mv TMP/jeter2.vcf.gz.tbi  TMP/jeter.vcf.gz.tbi
bcftools query -f '.' TMP/jeter.vcf.gz |wc -l 1>&2

#find frequent gnomad
bcftools query -f '%CHROM\t%POS0\t%END\\n'  TMP/jeter.vcf.gz > TMP/jeter.bed

if test ! -s TMP/jeter.bed
then
	echo 'chr22\t0\t1' > TMP/jeter.bed
fi

bcftools isec --regions-file TMP/jeter.bed -c none "${gnomad}" "TMP/jeter.vcf.gz" -O b --write-index -p TMP/isec1 -n =2 -w 1
bcftools view -O b  --write-index -i 'INFO/AF_nfe >= ${max_af }'  TMP/isec1/0000.bcf -o TMP/bad_gnomad.bcf
rm  TMP/isec1/0000.bcf  TMP/isec1/0000.bcf.csi

#remove frequent gnomad
bcftools isec -c none  "TMP/jeter.vcf.gz"  TMP/bad_gnomad.bcf  -O z -p TMP/isec2 -n ~10 -w 1
bcftools index -f -t TMP/isec2/0000.vcf.gz
mv TMP/isec2/0000.vcf.gz  TMP/jeter.vcf.gz
mv TMP/isec2/0000.vcf.gz.tbi  TMP/jeter.vcf.gz.tbi
bcftools query -f '.' TMP/jeter.vcf.gz |wc -l 1>&2



bcftools query -f '.' TMP/jeter.vcf.gz > TMP/flag.txt

if test -s TMP/flag.txt
then
	mv TMP/jeter.vcf.gz ${meta.sample}.${meta.type}.vcf.gz
	mv TMP/jeter.vcf.gz.tbi ${meta.sample}.${meta.type}.vcf.gz.tbi
fi
"""

stub:
    def prefix ="${meta.sample}.${meta.type}"
"""
touch versions.yml ${prefix}.vcf.gz ${prefix}.vcf.gz.tbi
"""
}



process CONCAT_SITES {
tag "${meta.id}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta),path("genotyped??.vcf.gz")
output:
	tuple val(meta),path("sites.vcf.gz"),path("sites.vcf.gz.tbi"),emit:vcf
script:
"""
set -x
mkdir -p TMP
find ./ -name "*.vcf.gz" | while read F
do
	bcftools view -G --header-only "\${F}" > TMP/vcf.header
	
	bcftools annotate -x 'ID,QUAL,FILTER,INFO'  "\${F}" |\\
		bcftools view -G --no-header >> TMP/vcf.body
	
done

sort  -T TMP TMP/vcf.body | uniq > TMP/vcf.body2

cat TMP/vcf.header TMP/vcf.body2 |\\
	bcftools sort --max-mem ${task.memory.giga}G -T TMP/x -O u '-' |\\
	bcftools norm  --multiallelics +both -O u '-' |\\
	bcftools view -m 2 -M 2 --type snps -O z -o TMP/jeter.vcf.gz
	
bcftools index -f -t TMP/jeter.vcf.gz


mv TMP/jeter.vcf.gz sites.vcf.gz
mv TMP/jeter.vcf.gz.tbi sites.vcf.gz.tbi
"""
stub:
    prefix = "sites"
"""
touch versions.yml ${prefix}.vcf.gz ${prefix}.vcf.gz.tbi
"""
}



process	GENOTYPE_CALL {
tag "${meta.id}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta4),path(sitesvcf),path(sitesvcftbi)
	tuple val(meta),path(bam),path(bai)
output:
   	tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
script:
   def jvm = task.ext.jvm?:"-XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -DGATK_STACKTRACE_ON_USER_EXCEPTION=true"
   def prefix = "${meta.id}.${meta.sample}.${bam.name}"
"""
mkdir -p TMP
set -x


gatk --java-options "${jvm}" HaplotypeCaller \\
	-L "${sitesvcf}" \\
	--alleles "${sitesvcf}" \\
	-R "${fasta}" \\
	-I "${bam}" \\
	-O "TMP/jeter.vcf.gz" \\
	 --output-mode EMIT_ALL_ACTIVE_SITES

bcftools annotate -x 'ID,INFO,FILTER,QUAL' -O z -o TMP/jeter2.vcf.gz TMP/jeter.vcf.gz
bcftools index -f -t TMP/jeter2.vcf.gz

mv TMP/jeter2.vcf.gz ${prefix}.vcf.gz
mv TMP/jeter2.vcf.gz.tbi ${prefix}.vcf.gz.tbi
"""
stub:
   def prefix = "${meta.id}.${meta.sample}.${bam.name}"
"""
touch versions.yml ${prefix}.vcf.gz ${prefix}.vcf.gz.tbi
"""
}


process BCFTOOLS_MERGE {
tag "${meta.id}"
label "process_short"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
memory "10G"
input:
	tuple val(meta),path("VCFS/*")
output:
	tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
script:
"""
mkdir -p TMP
find VCFS/ -name "*.vcf.gz" | sort -V > TMP/jeter.list
bcftools merge --force-samples --file-list TMP/jeter.list --missing-to-ref -O z -o merged.vcf.gz
bcftools index -f -t merged.vcf.gz
"""

stub:
   def prefix = "merged"
"""
touch versions.yml ${prefix}.vcf.gz ${prefix}.vcf.gz.tbi
"""
}


process CONCORDANCE {
tag "${meta.id} ${meta.sample}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta),path(truth_vcf),path(tbi)
	tuple val(meta2),path(merged),path(mergedidx)
output:
	tuple val(meta),path("*.txt"),emit:concordance
	path("versions.yml"),emit:versions
script:
    def prefix = "${meta.id}.${meta.sample}.concordance"
"""
mkdir -p TMP
bcftools gtcheck -u GT -g "${truth_vcf}" "${merged}" |\
	grep  '^DCv2' | cut -f2- |\\
	awk '{printf("%s\t%f\\n",\$0,(int(\$6)/(1.0*int(\$5))));}'  |\
	LC_ALL=C sort -T TMP -t '\t' -k7,7gr |\\
	awk '{printf("%s\t%s\\n",\$0,(NR==1?"best":"."));}'  > ${prefix}.txt

touch versions.yml
"""
stub:
    def prefix = "${meta.id}.${meta.sample}.concordance"
"""
touch versions.yml ${prefix}.txt
"""
}

process PLOT_CONCORDANCE {
label "process_single"
tag "${meta.id} ${meta.sample}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
        tuple val(meta),path(metrics)
		
output:
        tuple val(meta),path("*.pdf"),emit:pdf
		path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:"${meta.id}.${meta.sample}.concordances"
	def treshold = (task.ext.treshold?:0.7) as double
"""
cat << '__EOF__' | R --vanilla --no-save
table <- read.table("${metrics}",sep="\t",header=FALSE,stringsAsFactors=FALSE);
colnames(table) <-c("S1","S2","x1","x2","x3","GT_ALT","CONCORDANCE","x5")

maxc = max(table\$CONCORDANCE)

pdf(paste0(((1.0-maxc)*100.0),"${prefix}.pdf"));
plot(
        x = table\$GT_ALT,
        y = table\$CONCORDANCE,
        xlab="Number_of_matching_genotypes",
        ylab="Concordance",
        ylim=c(0,1.0),
        main="Concordances ${meta.id}",
        sub=""
        );


best <- table[table\$CONCORDANCE == maxc | table\$CONCORDANCE > ${treshold},]
text(x = best\$GT_ALT, y= best\$CONCORDANCE, labels = best\$S1 , pos = 2 , cex = 0.5 , col= "green") 

abline(h=${treshold},col="blue")
dev.off();
__EOF__

touch versions.yml
"""
stub:
	def prefix = task.ext.prefix?:"${meta.id}.${meta.sample}.concordances"
"""
touch versions.yml ${prefix}.pdf
"""
}



process MERGE_CONCORDANCES {
label "process_single"
tag "${meta.id}"
afterScript "rm -rf TMP"
input:
	tuple val(meta),path("DATA/??.txt")
output:
	tuple val(meta),path("*.tsv"),emit:table
	path("versions.yml"),emit:versions
script:
	def prefix= task.ext.prefix?:"${meta.id}.concordance"
"""
mkdir -p TMP
echo -e 'Query_Sample\tGenotyped_Sample\tDiscordance\tAverage_-log_P(HWE)\tNumber_of_sites_compared\tNumber_of_matching_genotypes\tfraction\tbest' > ${prefix}.tsv
find DATA -name "*.txt" -exec cat '{}' ';'  | LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k7,7gr  >> ${prefix}.tsv
touch versions.yml
"""
stub:
	def prefix= task.ext.prefix?:"${meta.id}.concordance"
"""
touch versions.yml ${prefix}.tsv
"""
}