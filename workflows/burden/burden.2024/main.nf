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
			map{it.containsKey("p_assoc")?it:it.plus(p_assoc:0.01)}


		conditions_ch.view{"CONDITION: $it"}

		X1 = xgene_ch.output.splitCsv(sep:'\t',header:true).
			combine(conditions_ch).
			filter{it[0].gene_name.equals("KCNH2")}.
			view{"$it"}

		eff=DOWNLOAD_SNPEFF_DB(genomeId)

		PER_GENE(genomeId, samplesheet,eff.output, vcf2bedcontig_ch,X1)
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
	LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n >> TMP/tmp2.bed




mv TMP/tmp2.bed genes.bed
"""
}


String getSnpeffDB(genomeId) {
	if(genomeId.equals("hs38me")) return "GRCh38.105";
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
tag "${gene.gene_name} ${condition.id}"
memory "3g"
input:
	val(genomeId)
	path(samplesheet)
	path(snpeffDir)
	path(vcf2bed)
	tuple val(gene),val(condition)
script:
	def genome = params.genomes[genomeId]
"""
module load snpeff/5.2
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
# min max alleles
#
bcftools view --min-alleles 2 --max-alleles ${condition.max_alleles} -O b -o TMP/jeter2.bcf TMP/jeter1.bcf
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

if ${!condition.gnomad.population.isEmpty() && condition.gnomad.AF>0}
then
	bcftools view TMP/jeter1.bcf |\\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar vcfgnomad \\
		--bufferSize 10000 \\
		--gnomad "${genome.gnomad_genome}" \\
		--fields "${condition.gnomad.population}"  \\
		--max-af '${condition.gnomad.AF}' |\
	bcftools view --apply-filters '.,PASS' -O b -o TMP/jeter2.bcf
        mv TMP/jeter2.bcf TMP/jeter1.bcf
fi

echo "${gene.gene_id}" > TMP/gene.id

bcftools view TMP/jeter1.bcf |\\
java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffiltergenes -g TMP/gene.id > TMP/jeter.vcf

java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfgenesplitter \\
		-E 'ANN/GeneId ANN/FeatureId VEP/GeneId VEP/Ensp VEP/Feature' \\
		--manifest TMP/jeter.mf \\
		-o \${PWD}/TMP \\
		TMP/jeter.vcf

cat << 'EOF' > TMP/jeter01.R
T1<-read.table("TMP/jeter.mf", header=TRUE, comment.char = "")
write.table(data.frame(splitter=T1\$splitter,key=T1\$key,path=T1\$path),quote=FALSE,sep="\\t",file=stdout(),col.names = FALSE,row.names = FALSE)
EOF

R --vanilla --slave < TMP/jeter01.R | while read SPLITTER KEY VCF
do

	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${HOME}/jvarkit.jar vcf2r \\
		--cases TMP/cases.txt \\
		--controls TMP/controls.txt "\${VCF}" > TMP/jeter.R

done


"""
}



process EXTRACT_EXONS {
memory "2g"
afterScript "rm -rf TMP"
input:
	val(genomeId)
	path(genesbed)
output:
	path("exons.bed"),emit:bed
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
	def gtf = genome.gtf
	def slop= params.genes_slop
"""
hostname 1>&2
${moduleLoad("htslib jvarkit bedtools")}
set -o pipefail
mkdir -p TMP

# remove header
awk -F '\t' '(\$1!="chrom")' '${genesbed}' > TMP/genes.bed

tabix --regions TMP/genes.bed  "${gtf}" |\
	java -jar \${JVARKIT_DIST}/gtf2bed.jar -R "${reference}" --columns "gtf.feature"|\
	awk -F '\t' '(\$4=="exon")' |\
	cut -f1,2,3 |\
	bedtools slop -i - -g "${reference}.fai" -b ${slop} |\
	LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n |\
	bedtools merge > TMP/tmp1.bed

mv TMP/tmp1.bed exons.bed

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">extract protein_coding genes from gtf</entry>
	<entry key="gtf">${gtf}</entry>
	<entry key="versions">${getVersionCmd("bedtools awk jvarkit/gtf2bed tabix")}</entry>
	<entry key="bed">${genesbed}</entry>
</properties>
EOF
"""
}



process RVTEST_PER_GENES {
tag "${L.collect{T->T.gene_name}.join(",")}"
afterScript "rm -rf TMP"
memory {task.attempt==1?'5G':'20G'}
errorStrategy 'retry'
maxRetries 3
input:
	val(meta)
	val(genomeId)
	path(pedigree)
	path(vcfbed)
	val(L)
errorStrategy "retry"

output:
	path("assoc.list"),emit:output
	path("version.xml"),emit:version
script:
	if(!params.containsKey("rvtests")) throw new IllegalArgumentException("params.rvtest");
	if(!params.rvtests.containsKey("arguments")) throw new IllegalArgumentException("params.rvtest.arguments missing");
	def rvtest_arguments = params.rvtests.arguments
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
"""
hostname 1>&2
${moduleLoad("rvtests bcftools jvarkit bedtools")}

mkdir -p TMP/ASSOC ASSOC

i=1


## https://github.com/zhanxw/rvtests/issues/80
cut -f 1 "${reference}.fai" > TMP/chroms.A.txt
sed 's/^chr//' TMP/chroms.A.txt > TMP/chroms.B.txt
paste TMP/chroms.A.txt TMP/chroms.B.txt > TMP/chroms.C.txt


# sort VCF bed
sort -t '\t' -T TMP -k1,1 -k2,2n  '${vcfbed}' > TMP/vcfs.bed


cat << EOF > TMP/input.tsv
${L.collect(T->T.chrom+"\t"+T.start+"\t"+T.end+"\t"+T.gene_name+"\t"+T.gene_id).join("\n")}
EOF

cat TMP/input.tsv | while IFS=\$'\\t' read CHROM START END GENE_NAME GENE_ID
do

# create a bed for this gene
echo "\${CHROM}\t\${START}\t\${END}" > TMP/gene.bed


# find VCF overlapping this gene
bedtools intersect -wa -u -a TMP/vcfs.bed -b TMP/gene.bed | cut -f4 | sort | uniq > TMP/vcfs.list


# vcf is empty (wgselect removed al the genes) . peek a random one
if test ! -s TMP/vcfs.list ; then
	head -n 1 TMP/vcfs.bed | cut -f 4 > TMP/vcfs.list
fi
test -s TMP/vcfs.list

# file containing the name of the gene
echo "\${GENE_ID}" > TMP/gene.id



bcftools concat --regions-file TMP/gene.bed --file-list TMP/vcfs.list -O u --allow-overlaps --remove-duplicates |\
	bcftools annotate -O v --rename-chrs TMP/chroms.C.txt |\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffiltergenes -g TMP/gene.id > TMP/jeter.vcf

	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfgenesplitter \
		-E 'ANN/GeneId ANN/FeatureId VEP/GeneId VEP/Ensp VEP/Feature' \
		--manifest TMP/jeter.mf \
		-o \${PWD}/TMP \
		TMP/jeter.vcf

	rm TMP/jeter.vcf

	awk -F '\t' '/^[^#]/ {S=\$4;gsub(/[\\/]+/,"_",S);G=\$5;gsub(/[\\/_\\-]+/,"_",G);K=\$6;gsub(/[\\/]+/,"_",K);printf("%s_%s_%s\t%s:%d-%d\t%s\\n",K,G,S,\$1,\$2,\$3,\$7);}' TMP/jeter.mf |\
	while read K RANGE V
	do
		bcftools index -t "\${V}"

		# build setFile
		echo "\$K\t\${RANGE}" > TMP/variants.setfile

		rvtest  --noweb \
        		--inVcf "\${V}" \
			--setFile TMP/variants.setfile \
	        	--pheno "${pedigree}" \
		        --out "TMP/ASSOC/part.\${GENE_ID}.0\${i}" \
			${rvtest_arguments} 2> TMP/last.rvtest.log

		rm "\${V}" "\${V}.tbi"

		i=\$((i+1))
	done

done

# collect each assoc by type and merge
find \${PWD}/TMP/ASSOC -type f -name "part.*assoc" | awk -F '/' '{N=split(\$NF,a,/\\./); print a[N-1];}' | sort -T TMP | uniq | while read SUFFIX
do

	find \${PWD}/TMP/ASSOC -type f -name "part.*\${SUFFIX}.assoc" > TMP/assoc.list
	test -s TMP/assoc.list

	xargs -a TMP/assoc.list cat | sed 's/^Range/#Range/' | LC_ALL=C sort -T TMP | uniq | sed 's/^#Range/Range/' > "ASSOC/concat.\${SUFFIX}.assoc"
done


find \${PWD}/ASSOC -type f -name "concat.*.assoc" >  assoc.list


cat << EOF > version.xml
<properties id="${task.process}">
  <entry key="name">${task.process}</entry>
  <entry key="description">VCF is split by transcript / gene, and then rvtest is applied.</entry>
  <entry key="name"Pedigree">${pedigree}</entry>
  <entry key="rvtest.path">\$(which rvtest)</entry>
  <entry key="rvtest.version">\$(rvtest  --version 2> /dev/null  | head -n1)</entry>
</properties>
EOF
"""
stub:
"""
touch assoc.list
echo "<properties/>" > version.xml
"""
}


