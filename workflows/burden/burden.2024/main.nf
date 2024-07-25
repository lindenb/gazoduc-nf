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

include {runOnComplete;} from '../../../modules/utils/functions.nf'


// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run my_pipeline --input input_file.csv")
   exit 0
}
// validate parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)



workflow {
		ch_input = Channel.fromList(samplesheetToList(params.samplesheet, "assets/schema_samplesheet.json"))
		BURDEN_CODING(Channel.fromPath(params.vcfs), file(params.samplesheet), Channel.fromPath(params.conditions), file(params.bed) )
		}

workflow BURDEN_CODING {
	take:
		vcfs
		samplesheet
		conditions
		bed
	main:

		ch1 = vcfs.
			splitText().
			map{it.trim()}.
			map{
			if(!(it.endsWith(".bcf") || it.endsWith(".vcf.gz")) ) throw new IllegalArgumentException("bad vcf extension for ${it}"); 
			if(it.endsWith(".bcf")){
				return [it,it+".csi"];
				}
			else if(it.vcf.endsWith(".vcf.gz")){
				return [it,it+".tbi"];
				}
			else	{
				throw new IllegalArgumentException("bad vcf extension ${it.vcf}");
				}
			}.
			map{[file(it[0]),file(it[1])]}

		genome_ch = Channel.value([file(params.fasta),file(params.fasta+".fai"),file(""+file(params.fasta).getParent()+"/"+file(params.fasta).getBaseName()+".dict")])
		gtf_ch = Channel.value([file(params.gtf),file(params.gtf+".tbi")])

		

		
		vcf2bed_ch = BCFTOOLS_INDEX_S(ch1)
		ctg_start_end1  = vcf2bed_ch.output.splitCsv(sep:'\t',header:false).
			map{tuple(it[0][0],it[0][1],it[0][2],it[1],it[2])}

		ctg_start_end = ctg_start_end1

		vcf2bedcontig_ch = ctg_start_end.
			map{rec->rec.join("\t")}.
			collectFile(name: 'vcf2bed.bed', newLine: true)



		conditions_ch = conditions.splitJson().
			map{it.containsKey("description")?it:it.plus(description:it.id)}.
			map{it.containsKey("enabled")?it:it.plus(enabled:true)}.
			map{it.containsKey("f_missing")?it:it.plus(f_missing:-1.0)}.
			map{it.containsKey("minDP")?it:it.plus(minDP:5)}.
			map{it.containsKey("maxDP")?it:it.plus(maxDP:300)}.
			map{it.containsKey("lowGQ")?it:it.plus(lowGQ:60)}.
			map{it.containsKey("minGQsingleton")?it:it.plus(minGQsingleton:90)}.
			map{it.containsKey("minRatioAD")?it:it.plus(minRatioAD:0.2)}.
			map{it.containsKey("so_acn")?it:it.plus(so_acn:"")}.
			map{it.containsKey("gnomad")?it:it.plus(gnomad:[:])}.
			map{
				if(it.gnomad.containsKey("population")) return it;
				it.gnomad.population="";
				return it;
			}.
			map{
				if(it.gnomad.containsKey("AF")) return it;
				it.gnomad.AF=0.0;
				return it;
			}.
			map{it.containsKey("max_alleles")?it:it.plus(max_alleles:3)}.
			map{it.containsKey("maxAF")?it:it.plus(maxAF:0.1)}.
			map{it.containsKey("p_assoc")?it:it.plus(p_assoc:0.001)}.
			map{it.containsKey("cadd")?it:it.plus(cadd:[:])}.
			map{
				if(it.cadd.containsKey("phred")) return it;
				it.cadd.phred=-1.0;
				return it;
			}.
			map {
				if(!it.containsKey("ignore_HOM_VAR")) it=it.plus("ignore_HOM_VAR":false);
				if(!it.containsKey("phastCons")) it=it.plus("phastCons":-1);
				if(!it.containsKey("bed_is_data_source")) it=it.plus("bed_is_data_source":false);
				if(!it.containsKey("types")) it=it.plus("types":"");
				if(!it.containsKey("gene_biotype")) it=it.plus("gene_biotype":"");
				if(!it.containsKey("polyx")) it=it.plus("polyx":10);
				if(!it.containsKey("exons_only")) it=it.plus("exons_only":false)
				if(!it.containsKey("QD")) it=it.plus("QD":2);
				if(!it.containsKey("FS")) it=it.plus("FS":60);
				if(!it.containsKey("SOR")) it=it.plus("SOR":3);
				if(!it.containsKey("MQ")) it=it.plus("MQ":40);
				if(!it.containsKey("MQRankSum")) it=it.plus("MQRankSum":-12.0);
				if(!it.containsKey("ReadPosRankSum")) it=it.plus("ReadPosRankSum":-8.0);
				return it;
			}
		

		conditions_ch.view{"CONDITION: $it"}


		xgene_ch = EXTRACT_GENES(genome_ch,gtf_ch,bed,vcf2bedcontig_ch.combine(conditions_ch))
		
		bed_cond_vcf_ch = xgene_ch.output.
			flatMap{X->X[0].flatten().collect{Y->tuple(Y,X[1])}}.
			combine(vcf2bedcontig_ch)
		
		gene_cond_vcf_ch = Channel.empty()


		eff=DOWNLOAD_SNPEFF_DB()

		per_gene_ch = PER_BED(
			genome_ch, 
			gtf_ch,
			samplesheet,
			(params.phastCons.equals("NO_FILE")? file("NO_PHASTCONS"): file(params.phastCons))  ,
			(params.excludeBed.equals("NO_FILE")? file("NO_EXCLUDE"): file(params.excludeBed))  ,
			eff.output,
			bed_cond_vcf_ch
			)
		per_gene_ch.output.splitCsv(sep:'\t', header:true).
			view{"$it"}
		all_results_ch =  per_gene_ch.output.collectFile(name:"output.tsv")
	
		cleanup_ch = CLEANUP(all_results_ch)

		CONCAT_VCFS(per_gene_ch.vcfs.flatten().collect())

		plot_ch=PLOTIT( 
			Channel.of("fisher","skat","skato","saktadj","saktoadj").
				combine(conditions_ch).
				combine(cleanup_ch.output)
			)
		MULTIQC(plot_ch.collect())

	}

process BCFTOOLS_INDEX_S {
tag "${vcf.name}"
label 'process_low'
input:
	tuple path(vcf),path(idx)
output:
	tuple path("${vcf.getSimpleName()}.bed"),path(vcf),path(idx),emit:output
script:
"""
bcftools index -s "${vcf}" | awk -F '\t' '{printf("%s\\t0\\t%s\\n",\$1,\$2);}' > "${vcf.getSimpleName()}.bed"
"""
}


process EXTRACT_GENES {
tag "${vcf2bed} ${condition.id}"
memory "2g"
afterScript "rm -rf TMP"
input:
	tuple path(fasta),path(fai),path(dict)
	tuple path(gtf), path(gtf_tbi)
	path(userbed)
	tuple path(vcf2bed),val(condition)
output:
	tuple path("BEDS/*.bed"),val(condition),emit:output
script:
	def slop = 20
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP

cut -f1,2,3 '${vcf2bed}' | LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n > TMP/tmp1.bed

#
# FOR EXAMPLE EACH RECORD IN THE BED IS A PEAK OF ATAC-SEQ
#
if ${condition.bed_is_data_source == true}
then

	# user bed MUST be defined
	${!userbed.name.equals("NO_FILE")}

	${userbed.name.endsWith(".gz")?"gunzip -c":"cat"} ${userbed} |\\
	cut -f1,2,3 |\\
	LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n |\\
	bedtools intersect -u -wa -a - -b TMP/tmp1.bed |\\
	LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n -k3,3n|\\
	uniq > TMP/tmp3.bed

else


tabix --regions "TMP/tmp1.bed" "${gtf}" |\\
	java -jar \${JVARKIT_DIST}/gtf2bed.jar -R "${fasta}" --columns "gtf.feature,gene_biotype" |\\
	awk -F '\t' '(\$4=="gene" ${condition.gene_biotype.isEmpty()?"":"&& \$5==\"${condition.gene_biotype}\""})' |\\
	cut -f1,2,3 |\\
	bedtools slop -i -  -g "${fai}" -b ${slop} |\\
	LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n |\\
	bedtools intersect -u -wa -a - -b TMP/tmp1.bed |\\
	LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n > TMP/tmp3.bed


if ${!userbed.name.equals("NO_FILE")}
then
	cut -f1,2,3 "${userbed}" |\
		LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n |\\
		bedtools merge > TMP/tmp4.bed

	bedtools intersect -u -wa -a TMP/tmp3.bed -b TMP/tmp4.bed |\\
		LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n >> TMP/tmp5.bed
	mv TMP/tmp5.bed TMP/tmp3.bed
fi

fi

mkdir -p BEDS

java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/bedcluster.jar \\
		-R "${fasta}" --out BEDS  --size "${params.group_by_cluster_size}" \\
		TMP/tmp3.bed
"""
}



process  DOWNLOAD_SNPEFF_DB {
memory "3G"
output:
	path("SNPEFF"),emit:output
script:
"""
module load snpeff/5.2
mkdir -p SNPEFFX TMP
snpEff -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  download -dataDir  "\${PWD}/SNPEFFX"  '${params.snpeff_db}'
test -s SNPEFFX/*/snpEffectPredictor.bin
mv SNPEFFX SNPEFF

"""
}

String countVariants(def f) {
	return "bcftools query -f '.\\n' \""+f+"\" | wc -l 1>&2";
	}


process PER_BED {
afterScript "rm -rf TMP"
tag "${condition.id} ${roi.name}"
memory "3g"
//array 100
errorStrategy "retry"
maxRetries 2
input:
	tuple path(fasta),path(fai),path(dict)
	tuple path(gtf), path(gtf_tbi)
	path(samplesheet)
	path(phastCons)
	path(excludeBed)
	path(snpeffDir)
	tuple path(roi),val(condition),path(vcf2bed)
output:
	path("results.bed"),emit:output
	path("VCFS/*"),optional:true,emit:vcfs
when:
	condition.enabled==true
script:
	def slop=20
"""
hostname 1>&2
module load snpeff/5.2 R/3.6.0-dev
mkdir -p TMP VCFS
set -x

# compile Minikit
cat "${moduleDir}/Minikit.java" > TMP/Minikit.java
javac -d TMP -cp \${JVARKIT_DIST}/jvarkit.jar -sourcepath TMP  TMP/Minikit.java


cat << 'EOF' > TMP/jeter.R
T1 <- read.csv("${samplesheet}",header=TRUE)
T2 <- T1[T1\$status=='case',]\$sample
head(T2)
write.table(T2,"TMP/cases.txt",sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)
T2 <- T1[T1\$status=='control',]\$sample
head(T2)
write.table(T2,"TMP/controls.txt",sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)
EOF


R --no-save --vanilla < TMP/jeter.R

sort TMP/cases.txt | uniq > TMP/jeter.txt
mv TMP/jeter.txt TMP/cases.txt
sort TMP/controls.txt | uniq > TMP/jeter.txt
mv TMP/jeter.txt TMP/controls.txt


LC_ALL=C sort -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n "${roi}" > TMP/roi.bed


if ${condition.exons_only==true}
then
	tabix --regions TMP/roi.bed "${gtf}" |\\
        	java -jar \${JVARKIT_DIST}/gtf2bed.jar -R "${fasta}" --columns "gtf.feature" |\\
		awk -F '\t' '(\$4=="exon")' |\\
		cut -f1-3 |\\
	        bedtools slop -i -  -g "${fai}" -b ${slop} |\\
        	LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n |\\
		bedtools merge > TMP/exons.bed
	
	if ! test -s TMP/exons.bed
	then
		echo -e 'chr22\t0\t1' > TMP/exons.bed
	fi

	bedtools intersect -a TMP/roi.bed -b TMP/exons.bed |\\
		LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n |\\
		bedtools merge > TMP/jeter.bed
	
	mv -v TMP/jeter.bed TMP/roi.bed

	if ! test -s TMP/roi.bed
	then
		echo -e 'chr22\t0\t1' > TMP/roi.bed
	fi


fi



bedtools intersect -u -a ${vcf2bed} -b TMP/roi.bed | cut -f 4 | sort | uniq > TMP/vcf.list

if ! test -s TMP/vcf.list
then
	head -n 1  ${vcf2bed} | cut -f 4 > TMP/vcf.list
fi

bcftools concat -a --regions-file TMP/roi.bed --file-list TMP/vcf.list -O b -o TMP/jeter1.bcf

${countVariants("TMP/jeter1.bcf")}

#
# BLACKLIST BED
#
if ${!excludeBed.name.equals("NO_EXCLUDE")}
then

	bcftools view --targets-file "^${excludeBed}" --targets-overlap 2  -O b -o TMP/jeter2.bcf TMP/jeter1.bcf
	mv TMP/jeter2.bcf TMP/jeter1.bcf
	${countVariants("TMP/jeter1.bcf")}

fi

bcftools query -l TMP/jeter1.bcf | sort |\
	comm -12 - TMP/cases.txt > TMP/jeter.txt
mv TMP/jeter.txt TMP/cases.txt

test -s TMP/cases.txt

bcftools query -l TMP/jeter1.bcf | sort |\
	comm -12 - TMP/controls.txt > TMP/jeter.txt
mv TMP/jeter.txt TMP/controls.txt

test -s TMP/controls.txt

comm -12 TMP/controls.txt TMP/cases.txt > TMP/jeter.txt
test ! -s TMP/jeter.txt

cat TMP/controls.txt TMP/cases.txt > TMP/all.samples.txt


bcftools view --trim-unseen-allele --trim-alt-alleles --samples-file TMP/all.samples.txt -O b -o TMP/jeter2.bcf TMP/jeter1.bcf
mv TMP/jeter2.bcf TMP/jeter1.bcf
${countVariants("TMP/jeter1.bcf")}


## too many alleles
bcftools view -m 2 -M ${condition.max_alleles} -O b -o TMP/jeter2.bcf TMP/jeter1.bcf
mv TMP/jeter2.bcf TMP/jeter1.bcf
${countVariants("TMP/jeter1.bcf")}

#
# normalize
#
bcftools norm -f ${params.fasta} --multiallelics -any -O u TMP/jeter1.bcf |\\
	bcftools view  -i 'ALT!="*" && AC[*]>0 && INFO/AF <= ${condition.maxAF}' -O b -o TMP/jeter2.bcf
mv TMP/jeter2.bcf TMP/jeter1.bcf
${countVariants("TMP/jeter1.bcf")}


#
# variant type
#
if ${!condition.types.isEmpty()}
then
	bcftools view --types "${condition.types}" -O b -o TMP/jeter2.bcf TMP/jeter1.bcf
	mv TMP/jeter2.bcf TMP/jeter1.bcf
	${countVariants("TMP/jeter1.bcf")}
fi



#
# bcftools contrast
#
bcftools +contrast -0 TMP/controls.txt -1 TMP/cases.txt -O u TMP/jeter1.bcf |\\
	bcftools view -e 'INFO/PASSOC < ${condition.p_assoc}' -O b -o TMP/jeter2.bcf
mv -v TMP/jeter2.bcf TMP/jeter1.bcf
${countVariants("TMP/jeter1.bcf")}



## polyx
if ${ condition.polyx  > 1 } ; then

	bcftools view TMP/jeter1.bcf |\\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar vcfpolyx -R "${fasta}" --tag POLYX -n "${condition.polyx}"  |\
                bcftools view -e 'INFO/POLYX > ${condition.polyx}' -O b -o  TMP/jeter2.bcf
	
	mv -v TMP/jeter2.bcf TMP/jeter1.bcf
	${countVariants("TMP/jeter1.bcf")}

fi


if ${!condition.so_acn.isEmpty()}
then
	bcftools view TMP/jeter1.bcf |\\
	snpEff -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  eff -dataDir "\${PWD}/${snpeffDir}" \\
		 -nodownload -noNextProt -noMotif -noInteraction -noLog -noStats -chr chr -i vcf -o vcf '${params.snpeff_db}' |\\
	java -jar \${JVARKIT_DIST}/jvarkit.jar vcffilterso \\
		--remove-attribute  --rmnoatt \\
		 -A '${condition.so_acn}' |\\
	bcftools view -O b -o TMP/jeter2.bcf
	mv TMP/jeter2.bcf TMP/jeter1.bcf
	${countVariants("TMP/jeter1.bcf")}

fi


if ${!phastCons.name.equals("NO_PHASTCONS") && condition.phastCons >= 0}
then
	bcftools view TMP/jeter1.bcf |\\
		java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar vcfbigwig -B '${phastCons}' --tag "PHASTCONS" | \\
		 bcftools view -e 'INFO/PHASTCONS < ${condition.phastCons}'  -O b -o TMP/jeter2.bcf

	mv TMP/jeter2.bcf TMP/jeter1.bcf
	${countVariants("TMP/jeter1.bcf")}

fi

#
# CADD
#
if ${params.containsKey("cadd") && (condition.cadd.phred as int) >0}
then
	bcftools view TMP/jeter1.bcf |\\
		java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar vcfcadd --tabix '${params.cadd}' | \\
		bcftools view -e 'INFO/CADD_PHRED < ${condition.cadd.phred}' -O b -o TMP/jeter2.bcf
	mv -v TMP/jeter2.bcf TMP/jeter1.bcf
	${countVariants("TMP/jeter1.bcf")}

fi


if ${!condition.gnomad.population.isEmpty() && condition.gnomad.AF>0}
then
	bcftools view TMP/jeter1.bcf |\\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar vcfgnomad \\
		--bufferSize 10000 \\
		--gnomad "${params.gnomad}" \\
		--fields "${condition.gnomad.population}"  \\
		--max-af '${condition.gnomad.AF}' |\\
	bcftools view --apply-filters '.,PASS' -O b -o TMP/jeter2.bcf
        mv TMP/jeter2.bcf TMP/jeter1.bcf
	${countVariants("TMP/jeter1.bcf")}

fi


cat << 'EOF' >> TMP/jeter02.R

v1 <- applyCAST(genotypes,population,variants)
v2 <- applySKAT(genotypes,population,variants)
v3 <- applySKATO(genotypes,population,variants)
v4 <- applySKATAdjusted(genotypes,population,variants)
v5 <- applySKATOAdjusted(genotypes,population,variants)

cat(
	c(
		contig,
		gene.start,
		gene.end,
		"${condition.id}",
		splitter,
		gene.key,
		gene.name,
		n.variants,
		v1\$p.value,
		v2\$p.value,
		v3\$p.value,
		v4\$p.value,
		v5\$p.value,
		ODD_RATIO,
		CASES_ALT_COUNT,
		CASES_REF_COUNT,
		CTRLS_ALT_COUNT,
		CTRLS_REF_COUNT,
		variants.str,
		CASES_NAMES,
		CTRLS_NAMES
		),
	file="TMP/results.bed",
	sep = "\t",
	append = TRUE
)

cat("\\n",
   file="TMP/results.bed",
   sep="",
   append = TRUE
   )


EOF


#
# output header
#
cat << EOF | paste -sd '\t' > TMP/results.bed
contig
start
end
condition_id
splitter
key
gene_name
n_variants
fisher
skat
skato
saktadj
saktoadj
ODD_RATIO
CASE_ALT
CASE_REF
CTRL_ALT
CTRL_REF
variants
CASES_NAMES
CTRLS_NAMES
EOF


bcftools view -O v TMP/jeter1.bcf |\\
java  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -cp \${JVARKIT_DIST}/jvarkit.jar:TMP  Minikit \\
		--vcf-out VCFS \\
		${condition.bed_is_data_source?"--bed ${roi}":""} \\
		--QD "${condition.QD}" \\
		--FS "${condition.FS}" \\
		--SOR "${condition.SOR}" \\
		--MQ "${condition.MQ}" \\
		--MQRankSum "${condition.MQRankSum}" \\
		--ReadPosRankSum "${condition.ReadPosRankSum}" \\
		--minDP ${condition.minDP} \\
		--maxDP ${condition.maxDP} \\
		--minGQsingleton ${condition.minGQsingleton} \\
		--minRatioAD ${condition.minRatioAD} \\
		--f_missing  ${condition.f_missing} \\
		--lowGQ ${condition.lowGQ} \\
		${condition.ignore_HOM_VAR?"--ignore_HOM_VAR":""} \\
		--gtf "${gtf}" \\
		--header "${moduleDir}/karaka01.R" \\
		--body TMP/jeter02.R \\
		--cases TMP/cases.txt \\
		--controls TMP/controls.txt | \\
	R --vanilla --no-save --slave

mv -v TMP/results.bed ./

"""
}

process CLEANUP {
input:
	path(tsv)
output:
	path("cleanup.bed"),emit:output
script:
"""
awk -F '\t' '(NR==1 || \$1!="contig")' '${tsv}' > jeter.tsv
mv jeter.tsv cleanup.bed
"""
}

process CONCAT_VCFS {
input:
	path("VCFS/*")
output:
	tuple path("${params.prefix}concat.bcf"), path("${params.prefix}concat.bcf.csi"),emit:output
script:
"""
module load bcftools
bcftools concat -a -O b "${params.prefix}concat.bcf" VCFS/*.vcf.gz
bcftools index -f "${params.prefix}concat.bcf"
"""
}

process PLOTIT {
tag "${assoc} ${condition.id}"
input:
	tuple val(assoc),val(condition),path(results)
output:
	path("mqc.*"),emit:output
script:
	def head = 20
	def subtitle ="${condition.id}/${assoc}"
	def assoc_desc= assoc; //"${testDescription(assoc)}"
	def associd=  "${assoc}_${condition.id}"
"""
hostname 1>&2
module load R/3.6.0-dev
mkdir -p TMP


cat << '__EOF__' > jeter.R

T1<-read.table("${results}",header=TRUE,sep="\t",stringsAsFactors=FALSE)
head(T1)


T1<-T1[T1\$condition_id=="${condition.id}",]
head(T1)



T1<-T1[!is.na(as.numeric(as.character(T1\$${assoc}))),]
head(T1)


T2<-head(T1[order(as.numeric(T1\$${assoc})),],n=${head})
write.table(T2,file="TMP/jeter.tsv",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)


T1<-T1[,c("contig","start","key","${assoc}")]
head(T1)


T1\$contig = as.numeric(sub("Y","24",sub("X","23",sub("chr","",T1\$contig))))
head(T1)

T1<-T1[!is.na(T1\$contig),]
head(T1)





colnames(T1)<-c("CHR","BP","SNP","P")

library("qqman",lib.loc="/LAB-DATA/BiRD/users/lindenbaum-p/R")


if(nrow(T1)>0) {

Sys.setenv("DISPLAY"=":0.0")

png("mqc.${associd}.manhattan.png", type="cairo")
manhattan(T1,main="${associd}",sub="${subtitle}");
dev.off()

png("mqc.${associd}.qqplot.png", type="cairo")
qq(T1\$P,main="${associd}",sub="${subtitle}");
dev.off()


}
__EOF__

R --vanilla < jeter.R || true


cat << EOF > mqc.${associd}.yaml
custom_data:
  ${associd}_manhattan:
    parent_id: ${associd}_section
    parent_name: "${associd}"
    parent_description: "${assoc_desc} <pre>${condition}</pre>"
    section_name: "${associd} Manhattan"
    description: "${associd} Manhattan plot"
  ${associd}_qqplot:
    parent_id: ${associd}_section
    parent_name: "${associd}"
    parent_description: "${assoc_desc} <pre>${condition}</pre>"
    section_name: "${associd} QQPlot"
    description: "${associd} QQPlot"
sp:
  ${associd}_manhattan:
    fn: "mqc.${associd}.manhattan.png"
  ${associd}_qqplot:
    fn: "mqc.${associd}.qqplot.png"
ignore_images: false
EOF

cat << EOF > "mqc.${associd}.table_mqc.html"
<!--
parent_id: ${associd}_section
parent_name: "${associd}"
parent_description: "${associd}  <pre>${condition}</pre>"
id: '${associd}_table'
section_name: '${associd} table'
description: '${head} first lines.'
-->
EOF

awk -F '\t' 'BEGIN{printf("<table class=\\"table\\">");} (NR==1) {printf("<thead><caption>${associd}</caption><tr>");for(i=1;i<=NF;i++) {printf("<th>%s</th>",\$i);}printf("</tr></thead><tbody>"); next;} {printf("<tr>");for(i=1;i<=NF;i++) printf("<td>%s</td>",\$i); printf("</tr>");} END{if(NR>1) printf("</tbody>");printf("</table><br/>\\n");}' TMP/jeter.tsv >> "mqc.${associd}.table_mqc.html"

"""
}

process MULTIQC {
executor "local"
input:
	path(mqcs)
script:
	def prefix = params.prefix
	def title="";
	def comment="";
"""

		hostname 1>&2
		module load multiqc
		mkdir -p TMP

# do NOT use -type f, those are symlinks
find .  -name "mqc.*" | grep -v '\\.yaml\$' > TMP/jeter.list


		mkdir -p "${prefix}multiqc"

		export LC_ALL=en_US.utf8
		multiqc  --filename  "${prefix}multiqc_report.html" --no-ansi \
			${title}  \\
			${comment}  \\
			--force \\
			`find . -name "*.yaml" | awk '{printf(" --config  %s ",\$0);}' ` \\
			--outdir "${prefix}multiqc" \\
			--file-list TMP/jeter.list
		
		zip -9 -r "${prefix}multiqc.zip" "${prefix}multiqc"
"""
}
