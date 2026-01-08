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
include { GHOSTSCRIPT_MERGE             } from '../../../modules/gs/merge'
include { BCFTOOLS_QUERY_SAMPLES        } from '../../../modules/bcftools/query.samples'
workflow CONCORDANCE_BAMS {
    take:
        metadata
        fasta
        fai
        dict
        truth_sites_vcfs // meta,vcf,idx
        bams
    main:
        
        versions = Channel.empty()
        multiqc = Channel.empty()
    
        MAKE_SITE_ONLY_VCF(
            fasta,
            fai,
            truth_sites_vcfs
                .map{_meta,vcf,_tbi->vcf}
                .collect()
                .map{files->[[id:metadata.id],files.sort()]}
            )
        versions = versions.mix(MAKE_SITE_ONLY_VCF.out.versions)


        GENOTYPE_CALL(
            fasta,
            fai,
            dict,
            MAKE_SITE_ONLY_VCF.out.vcf,
            bams.filter{meta,_bam,_bai->meta.type!="truth"}
            )
        versions = versions.mix(GENOTYPE_CALL.out.versions)

        BCFTOOLS_MERGE(
            fasta,
            fai,
            GENOTYPE_CALL.out.vcf
                .map{_meta,vcf,tbi->[vcf,tbi]}
                .flatMap()
                .collect()
                .map{[[id: (metadata.id?:"concordance")],it.sort()]}
            )
        versions = versions.mix(BCFTOOLS_MERGE.out.versions)

        BCFTOOLS_QUERY_SAMPLES( truth_sites_vcfs)
        versions = versions.mix(BCFTOOLS_QUERY_SAMPLES.out.versions)
        samples_ch = BCFTOOLS_QUERY_SAMPLES.out.samples.splitText()
            .map{meta,sn->[meta,sn.trim()]}
        
    // .map{meta,vcf,tbi,sn->[meta.plus(id:sn),vcf,id]}
        CONCORDANCE(
            truth_sites_vcfs.combine(samples_ch)
                .filter{meta1,_vcf,_tbi,meta2,_sn-> meta1==meta2}
                .map{meta1,vcf,tbi,_meta2,sn->[meta1.plus(id:sn),vcf,tbi]},
            BCFTOOLS_MERGE.out.vcf
            )
        versions = versions.mix(CONCORDANCE.out.versions)

        PLOT_CONCORDANCE(CONCORDANCE.out.concordance)
        versions = versions.mix(PLOT_CONCORDANCE.out.versions)


        GHOSTSCRIPT_MERGE(
            PLOT_CONCORDANCE.out.pdf
                .map{_meta,f->f}
                //.map{[it.name,it]} //PREVENT FILE NAME COLLISTION
                //.groupTuple()
                //.map{f,L->L[0]}
                .collect()
                .map{files->[ [id:metadata.id],files.sort()]}
            )
        versions = versions.mix(GHOSTSCRIPT_MERGE.out.versions)

        MERGE_CONCORDANCES( 
            CONCORDANCE.out.concordance
                .map{_meta,f->f}
                .collect()
                .map{files->[  [id:(metadata.id?:"concordance")] ,files.sort()]}
            )
        versions = versions.mix(MERGE_CONCORDANCES.out.versions)
    emit:
        versions
        multiqc
    }


process MAKE_SITE_ONLY_VCF {
tag "${meta.id}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta),path("genotyped??.vcf.gz")
output:
	tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"${meta.id}.sites"
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


mv TMP/jeter.vcf.gz "${prefix}.vcf.gz"
mv TMP/jeter.vcf.gz.tbi "${prefix}.vcf.gz.tbi"
touch versions.yml
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
    path("versions.yml"),emit:versions
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
touch versions.yml
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
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
	tuple val(meta),path("VCFS/*")
output:
	tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
    path("versions.yml"),emit:versions
script:
    def args = task.ext.args?:"--missing-to-ref"
"""
mkdir -p TMP
find VCFS/ -name "*.vcf.gz" | sort -V > TMP/jeter.list
bcftools merge --force-samples --file-list TMP/jeter.list ${args} -O u |\\
    bcftools norm --rm-dup all --multiallelics +both --fasta-ref ${fasta} -O u |\\
    bcftools view -c 1 -m 2 -M 2 -O z -o merged.vcf.gz
bcftools index -f -t merged.vcf.gz

touch versions.yml
"""

stub:
   def prefix = "merged"
"""
touch versions.yml ${prefix}.vcf.gz ${prefix}.vcf.gz.tbi
"""
}


process CONCORDANCE {
tag "${meta.id} ${meta1.id}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(truth_vcf),path(tbi)
	tuple val(meta ),path(merged),path(mergedidx)
output:
	tuple val(meta1),path("*.txt"),emit:concordance
	path("versions.yml"),emit:versions
script:
    def sample = "${meta1.id}"
    def prefix = "${meta.id}.${meta1.id}.concordance"
"""
mkdir -p TMP

bcftools view --samples "${sample}"  -O u "${truth_vcf}" |\\
    bcftools view -c 1 -O u |\\
    bcftools +setGT  -O b  -o "TMP/jeter.bcf"  -  --  -t a -n u

bcftools index -f TMP/jeter.bcf

bcftools gtcheck -u GT -g  TMP/jeter.bcf "${merged}" |\
	grep  '^DCv2' | cut -f2- |\\
	awk  -F '\t' '{F=-1; if(int(\$5) > 0 ) F=(int(\$6)/(1.0*int(\$5)));  printf("%s\t%f\\n",\$0,F);}'  |\
	LC_ALL=C sort -S '${task.memory.kilo}' -T TMP -t '\t' -k2,2 -k7,7gr |\\
    awk -F '\t' '{printf("%s\t%s\\n",\$0,(NR==1 && int(\$5) ?"best":"."));}'  > ${prefix}.txt

touch versions.yml
"""
stub:
    def prefix = "${meta.id}.${meta1.id}.concordance"
"""
touch versions.yml ${prefix}.txt
"""
}

process PLOT_CONCORDANCE {
label "process_single"
tag "${meta.id} "
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
        tuple val(meta),path(metrics)
		
output:
        tuple val(meta),path("*.pdf"),emit:pdf
		path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:"${meta.id}.concordances"
	def treshold = (task.ext.treshold?:0.7) as double
"""
cat << '__EOF__' | R --vanilla --no-save
table <- read.table("${metrics}",sep="\t",header=FALSE,stringsAsFactors=FALSE);
colnames(table) <-c("S1","S2","x3","x4","x5","GT_ALT","CONCORDANCE","x5")

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
