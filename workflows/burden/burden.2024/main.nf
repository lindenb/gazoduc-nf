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


// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)


workflow {
		ch_input = Channel.fromList(samplesheetToList(params.samplesheet, "assets/schema_samplesheet.json"))
		BURDEN_CODING(params.genomeId, Channel.fromPath(params.datasheet), file(params.samplesheet), Channel.fromPath(params.conditions), file(params.bed) )
		}

workflow BURDEN_CODING {
	take:
		genomeId
		datasheet
		samplesheet
		conditions
		bed
	main:

		ch1 = datasheet.
			splitCsv(sep:',',header:true).
			/** check VCF suffix */
			map{
			if(!it.containsKey("vcf")) throw new IllegalArgumentException("vcf key missing in datasheet");
			if(!(it.vcf.endsWith(".bcf") || it.vcf.endsWith(".vcf.gz")) ) throw new IllegalArgumentException("bad vcf extension for ${it.vcf}"); 
			return it;
			}.map{
			if(it.containsKey("idx") && !(it.idx.isEmpty() || it.idx.equals("."))) {
				if(!(it.idx.endsWith(".csi") || it.vcf.endsWith(".tbi"))) throw new IllegalArgumentException("bad vcf index extension for ${it.vcf}"); 
				return it;
				}
			else if(it.vcf.endsWith(".bcf")){
				return it.plus(idx: it.vcf+".csi");
				}
			else if(it.vcf.endsWith(".vcf.gz")){
				return it.plus(idx: it.vcf+".tbi");
				}
			else	{
				throw new IllegalArgumentException("bad vcf extension ${it.vcf}");
				}
			}.
			branch{
				with_coord: it.containsKey("chrom") && it.containsKey("start") && it.containsKey("end")
				no_coord : true // A fallback condition
				}

		

		
		vcf2bed_ch = BCFTOOLS_INDEX_S(ch1.no_coord.map{tuple(file(it.vcf),file(it.idx))})
		ctg_start_end1  = vcf2bed_ch.output.splitCsv(sep:'\t',header:false).
			map{tuple(it[0][0],it[0][1],it[0][2],it[1],it[2])}

		ctg_start_end2 = ch1.with_coord.map{tuple(it.chrom,it.start,it.end,file(it.vcf),file(it.idx))}
		

		ctg_start_end = ctg_start_end1.mix(ctg_start_end2)

		vcf2bedcontig_ch = ctg_start_end.
			map{rec->rec.join("\t")}.
			collectFile(name: 'vcf2bed.bed', newLine: true)


		xgene_ch = EXTRACT_GENES(genomeId,vcf2bedcontig_ch,bed)

		conditions_ch = conditions.splitJson().
			map{it.containsKey("description")?it:it.plus(description:it.id)}.
			map{it.containsKey("gene_biotype_regex")?it:it.plus(gene_biotype_regex:".*")}.
			map{it.containsKey("enabled")?it:it.plus(enabled:true)}.
			map{it.containsKey("f_missing")?it:it.plus(f_missing:-1.0)}.
			map{it.containsKey("minDP")?it:it.plus(minDP:5)}.
			map{it.containsKey("maxDP")?it:it.plus(maxDP:300)}.
			map{it.containsKey("lowGQ")?it:it.plus(lowGQ:60)}.
			map{it.containsKey("minGQsingleton")?it:it.plus(minGQsingleton:90)}.
			map{it.containsKey("minRatioSingleton")?it:it.plus(minRatioSingleton:0.2)}.
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
			map{it.containsKey("p_assoc")?it:it.plus(p_assoc:0.01)}.
			map{it.containsKey("cadd")?it:it.plus(cadd:[:])}.
			map{
				if(it.cadd.containsKey("phred")) return it;
				it.cadd.phred=-1.0;
				return it;
			}.
			map {
				if(!it.containsKey("QD")) it=it.plus("QD":2);
				if(!it.containsKey("FS")) it=it.plus("FS":60);
				if(!it.containsKey("SOR")) it=it.plus("SOR":3);
				if(!it.containsKey("MQ")) it=it.plus("MQ":40);
				if(!it.containsKey("MQRankSum")) it=it.plus("MQRankSum":-12.0);
				if(!it.containsKey("ReadPosRankSum")) it=it.plus("ReadPosRankSum":-8.0);
				return it;
			}
		

		conditions_ch.view{"CONDITION: $it"}

		X1 = xgene_ch.output.splitCsv(sep:'\t',header:true).
			combine(conditions_ch).
			combine(vcf2bedcontig_ch)



		eff=DOWNLOAD_SNPEFF_DB(genomeId)

		per_gene_ch = PER_GENE(genomeId, samplesheet,eff.output,X1)
		per_gene_ch.output.splitCsv(sep:'\t', header:true).
			view{"$it"}
		all_results_ch =  per_gene_ch.output.collectFile(name:"output.tsv")
	
		PLOTIT( 
			Channel.of("fisher","skat","skato","saktadj","saktoadj").
				combine(conditions_ch).
				combine(all_results_ch)
			)

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
tag "${vcf2bed}"
memory "2g"
afterScript "rm -rf TMP"
input:
	val(genomeId)
	path(vcf2bed)
	path(userbed)
output:
	path("genes.bed"),emit:output
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
	def gtf = genome.gtf
	def slop= params.genes_slop
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP

cut -f1,2,3 '${vcf2bed}' | LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n > TMP/tmp1.bed


cat << EOF  | paste -sd'\t' > TMP/tmp2.bed
contig
start
end
feature
gene_id
gene_name
strand
biotype
EOF

tabix --regions "TMP/tmp1.bed" "${gtf}" |\
	java -jar \${JVARKIT_DIST}/gtf2bed.jar -R "${reference}" --columns "gtf.feature,gene_id,gene_name,gff.strand,gene_biotype" |\
	awk -F '\t' '(\$4=="gene" && \$5!=".")' |\
	bedtools slop -i -  -g "${reference}.fai" -b ${slop} |\
	LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n |\
	bedtools intersect -u -wa -a - -b TMP/tmp1.bed |\
	LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n >> TMP/tmp3.bed


if ${!userbed.name.equals("NO_FILE")}
then
	cut -f1,2,3 "${userbed}" |\
		LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n |\\
		bedtools merge > TMP/tmp4.bed

	bedtools intersect -u -wa -a TMP/tmp3.bed -b TMP/tmp4.bed |\\
		LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n >> TMP/tmp5.bed
	mv TMP/tmp5.bed TMP/tmp3.bed
fi


cat TMP/tmp2.bed TMP/tmp3.bed > genes.bed
"""
}


String getSnpeffDB(genomeId) {
	if(genomeId.equals("hs38me")) return "GRCh38.105";
	if(genomeId.equals("hs37d5")) return "GRCh37.75";
	return "UNDEFINED_SNPEFF_DB";
	}

process  DOWNLOAD_SNPEFF_DB {
memory "3G"
input:
	val(genomeId)
output:
	path("SNPEFF"),emit:output
script:
"""
module load snpeff/5.2
mkdir -p SNPEFFX TMP
snpEff -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  download -dataDir  "\${PWD}/SNPEFFX"  ${getSnpeffDB(genomeId)}
test -s SNPEFFX/*/snpEffectPredictor.bin
mv SNPEFFX SNPEFF

"""
}

process PER_GENE {
afterScript "rm -rf TMP"
tag "${gene.gene_name} ${condition.id}"
memory "3g"
array 100
errorStrategy "retry"
maxRetries 2
input:
	val(genomeId)
	path(samplesheet)
	path(snpeffDir)
	tuple val(gene),val(condition),path(vcf2bed)
output:
	path("results.txt"),emit:output
when:
	condition.enabled==true && gene.biotype.matches(condition.gene_biotype_regex)
script:
	def genome = params.genomes[genomeId]
"""
module load snpeff/5.2 R/3.6.0-dev
mkdir -p TMP

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


echo "${gene.contig}\t${gene.start}\t${gene.end}" > TMP/gene.bed

bedtools intersect -u -a ${vcf2bed} -b TMP/gene.bed | cut -f 4 > TMP/vcf.list

if ! test -s TMP/vcf.list
then
	head -n 1  ${vcf2bed} | cut -f 4 > TMP/vcf.list
fi

bcftools concat -a --regions-file TMP/gene.bed --file-list TMP/vcf.list -O b -o TMP/jeter1.bcf

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


#
# Filter with jvarkit
#
cat << __EOF__ > TMP/jeter.code

/** min max alleles */
if(variant.getNAlleles()<2 || variant.getNAlleles()>${condition.max_alleles}) {
	return false;
	}


final double f_missing = ${condition.f_missing};
if(f_missing >= 0 ) {
	/** missing */
	final double  n_missing = variant.getGenotypes().stream().
		filter(G->G.isNoCall() || (G.hasDP() && G.getDP()==0)).
		count();

	if(n_missing/variant.getNSamples() > f_missing) return false;
	}

/** filters 
if(variant.isFiltered()) {
	return false;
}
*/

/** low DP */
final double dp= variant.getGenotypes().stream().
	filter(G->G.isCalled() && G.hasDP()).
	mapToInt(G->G.getDP()).average().orElse(${condition.minDP}); 

if(dp <${condition.minDP} || dp > ${condition.maxDP}) return false;

/** low GQ */

final long count_alt = variant.getGenotypes().stream().
        filter(g->g.isCalled() && !g.isHomRef()).
	count();

if(count_alt==0L) return false;

final double count_low_gq = variant.getGenotypes().stream().
	filter(g->g.isCalled() && !g.isHomRef() && g.hasGQ()).
	filter(g->g.getGQ()<${condition.lowGQ}).
	count();

if(count_low_gq/count_alt >= 0.25) {
	return false;
	}

// see https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
if(variant.hasAttribute("QD") && variant.getAttributeAsDouble("QD",1000) < ${condition.QD}) return false;
if(variant.hasAttribute("FS") && variant.getAttributeAsDouble("FS",0) > ${condition.FS}) return false;
if(variant.hasAttribute("SOR") && variant.getAttributeAsDouble("SOR",0) > ${condition.SOR}) return false;
if(variant.hasAttribute("MQ") && variant.getAttributeAsDouble("MQ",1000) < ${condition.MQ}) return false;
if(variant.hasAttribute("MQRankSum") && variant.getAttributeAsDouble("MQRankSum",1000) < ${condition.MQRankSum}) return false;
if(variant.hasAttribute("ReadPosRankSum") && variant.getAttributeAsDouble("ReadPosRankSum",1000) < ${condition.ReadPosRankSum}) return false;



Genotype singleton=null;
for(final Genotype g: variant.getGenotypes()) {
	if(g.isCalled() && !g.isHomRef()) {
		if(singleton!=null) return true;
		singleton=g;
		}
	}
//if(singleton!=null && singleton.isFiltered()) return false;

if(singleton!=null && singleton.isHet() && singleton.hasGQ() && singleton.getGQ()<${condition.minGQsingleton}) {
	return false; 
	}
if(singleton !=null && singleton.hasAD() && singleton.isHet() && singleton.getAD().length==2) 
	{
	int array[]=singleton.getAD();double r= array[1]/(double)(array[0]+array[1]);
	if(r< ${condition.minRatioSingleton} || r>(1.0 - ${condition.minRatioSingleton})) return false;
	}
return true;
__EOF__


bcftools view TMP/jeter1.bcf |\\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffilterjdk --nocode  -f TMP/jeter.code |\
	bcftools view  -O b -o TMP/jeter2.bcf
 mv -v TMP/jeter2.bcf TMP/jeter1.bcf


#
# normalize
#
bcftools norm -f ${genome.fasta} --multiallelics -any -O u TMP/jeter1.bcf |\\
	bcftools view  -i 'ALT!="*" && AC[*]>0 && INFO/AF <= ${condition.maxAF}' -O b -o TMP/jeter2.bcf
mv TMP/jeter2.bcf TMP/jeter1.bcf

#
# bcftools contrast
#
bcftools +contrast -0 TMP/controls.txt -1 TMP/cases.txt -O u TMP/jeter1.bcf |\\
	bcftools view -e 'INFO/PASSOC < ${condition.p_assoc}' -O b -o TMP/jeter2.bcf
mv -v TMP/jeter2.bcf TMP/jeter1.bcf





if ${!condition.so_acn.isEmpty()}
then
	bcftools view TMP/jeter1.bcf |\\
	snpEff -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  eff -dataDir "\${PWD}/${snpeffDir}" \\
		 -nodownload -noNextProt -noMotif -noInteraction -noLog -noStats -chr chr -i vcf -o vcf  ${getSnpeffDB(genomeId)} |\\
	java -jar \${JVARKIT_DIST}/jvarkit.jar vcffilterso \\
		--remove-attribute  --rmnoatt \\
		 -A '${condition.so_acn}' |\\
	bcftools view -O b -o TMP/jeter2.bcf
	mv TMP/jeter2.bcf TMP/jeter1.bcf
fi


#
# select only in the current gene
#
echo "${gene.gene_id}" > TMP/gene.id

bcftools view TMP/jeter1.bcf |\\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffiltergenes -g TMP/gene.id |\
	bcftools view  -O b -o TMP/jeter2.bcf
 mv -v TMP/jeter2.bcf TMP/jeter1.bcf

#
# CADD
#
if ${genome.containsKey("cadd_tabix") && condition.cadd.phred >0}
then
	bcftools view TMP/jeter1.bcf |\\
		java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar vcfcadd --tabix '${genome.cadd_tabix}' | \\
		bcftools view -e 'INFO/CADD_PHRED<${condition.cadd.phred}' -O b -o TMP/jeter2.bcf
	mv -v TMP/jeter2.bcf TMP/jeter1.bcf
fi


if ${!condition.gnomad.population.isEmpty() && condition.gnomad.AF>0}
then
	bcftools view TMP/jeter1.bcf |\\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar vcfgnomad \\
		--bufferSize 10000 \\
		--gnomad "${genome.gnomad_genome}" \\
		--fields "${condition.gnomad.population}"  \\
		--max-af '${condition.gnomad.AF}' |\\
	bcftools view --apply-filters '.,PASS' -O b -o TMP/jeter2.bcf
        mv TMP/jeter2.bcf TMP/jeter1.bcf
fi




bcftools view -O v TMP/jeter1.bcf |\\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfgenesplitter \\
		-E 'ANN/GeneId ANN/FeatureId VEP/GeneId VEP/Ensp VEP/Feature' \\
		--manifest TMP/jeter.mf \\
		-o \${PWD}/TMP

cat << 'EOF' > TMP/jeter01.R
T1<-read.table("TMP/jeter.mf", header=TRUE, comment.char = "")
write.table(data.frame(splitter=T1\$splitter,key=T1\$key,path=T1\$path),quote=FALSE,sep="\\t",file=stdout(),col.names = FALSE,row.names = FALSE)
EOF

echo -e '' > TMP/result.tsv

cat "${moduleDir}/karaka01.R" > TMP/jeter02.R

cat << 'EOF' >> TMP/jeter02.R

v1 <- applyCAST(genotypes,population,variants)
v2 <- applySKAT(genotypes,population,variants)
v3 <- applySKATO(genotypes,population,variants)
v4 <- applySKATAdjusted(genotypes,population,variants)
v5 <- applySKATOAdjusted(genotypes,population,variants)

cat(
	c(
		v1\$p.value,
		v2\$p.value,
		v3\$p.value,
		v4\$p.value,
		v5\$p.value,
		v1\$CASE_ALT,
		v1\$CASE_REF,
		v1\$CTRL_ALT,
		v1\$CTRL_REF
		),
	file="TMP/results.txt",
	sep = "\t",
	append = TRUE
)

EOF

#
# output header
#
cat << EOF | paste -sd '\t' > TMP/results.txt
condition_id
contig
start
end
gene_id
gene_name
gene_biotype
splitter
key
fisher
skat
skato
saktadj
saktoadj
CASE_ALT
CASE_REF
CTRL_ALT
CTRL_REF
EOF

R --vanilla --slave < TMP/jeter01.R | while read SPLITTER KEY VCF
do
	echo "VCF=\${VCF}" 1>&2

	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcf2r \\
		--variants \\
		--info AF \\
		--cases TMP/cases.txt \\
		--controls TMP/controls.txt "\${VCF}" > TMP/jeter.R

	cat TMP/jeter02.R >> TMP/jeter.R
	echo -n "${condition.id}\t${gene.contig}\t${gene.start}\t${gene.end}\t${gene.gene_id?:"."}\t${gene.gene_name?:"."}\t${gene.biotype?:"."}\t\${SPLITTER}\t\${KEY}\t" >> TMP/results.txt
	R --vanilla --no-save --slave < TMP/jeter.R

	# add EOL
	echo >> TMP/results.txt

done


mv -v TMP/results.txt ./

"""
}



process PLOTIT {
tag "${assoc} ${condition.id}"
input:
	tuple val(assoc),val(condition),path(results)
script:
	prefix = params.prefix?:""
	def head = 20
	def subtitle ="${condition.id}/${assoc}"
	def assoc_desc= assoc; //"${testDescription(assoc)}"
"""
hostname 1>&2
module load R/3.6.0-dev
mkdir -p TMP

# remove header
awk -F '\t' '(NR==1 || \$1!="condition_id")' '${results}' > TMP/jeter.tsv

wc -l TMP/jeter.tsv
head TMP/jeter.tsv

cat << '__EOF__' > jeter.R

T1<-read.table("TMP/jeter.tsv",header=TRUE,sep="\t",stringsAsFactors=FALSE)
head(T1)

T1<-T1[,c("contig","start","key","${assoc}")]
head(T1)

T1<-T1[!is.na(as.numeric(as.character(T1\$${assoc}))),]
head(T1)

colnames(T1)<-c("CHR","BP","SNP","P")

library("qqman",lib.loc="/LAB-DATA/BiRD/users/lindenbaum-p/R")


if(nrow(T1)>0) {

Sys.setenv("DISPLAY"=":0.0")

png("${prefix}${assoc}.manhattan.png")
manhattan(T1,main="${prefix}${assoc}",sub="${subtitle}");
dev.off()

png("${prefix}${assoc}.qqplot.png")
qq(T1\$P,main="${prefix}${assoc}",sub="${subtitle}");
dev.off()
}
__EOF__

R --vanilla < jeter.R || true


cat << EOF > multiqc_config.yaml
custom_data:
  ${assoc}_manhattan:
    parent_id: ${assoc}_section
    parent_name: "${assoc}"
    parent_description: "${assoc_desc}"
    section_name: "${assoc} Manhattan"
    description: "RVTEST ${assoc} Manhattan plot"
  ${assoc}_qqplot:
    parent_id: ${assoc}_section
    parent_name: "${assoc}"
    parent_description: "RVTEST ${assoc}"
    section_name: "${assoc} QQPlot"
    description: "RVTEST ${assoc} QQPlot"
sp:
  ${assoc}_manhattan:
    fn: "${prefix}${assoc}.manhattan.png"
  ${assoc}_qqplot:
    fn: "${prefix}${assoc}.qqplot.png"
ignore_images: false
EOF

cat << EOF > "${prefix}${assoc}.table_mqc.html"
<!--
parent_id: ${assoc}_section
parent_name: "${assoc}"
parent_description: "RVTEST ${assoc}"
id: '${assoc}_table'
section_name: '${assoc} table'
description: '${head} first lines.'
-->
<pre>
EOF

head -n ${head} TMP/jeter.tsv | column -t >> "${prefix}${assoc}.table_mqc.html"
echo "</pre>" >>  "${prefix}${assoc}.table_mqc.html"



find \${PWD} -type f -name "*.png" >> paths.txt
"""
}

