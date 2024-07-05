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

/***


Original code from Floriane Simonet

*/


nextflow.enable.dsl=2


/**

step_id and step_name are used to get a unique id for multiqc when using this submodule multiple times

*/
if(!params.containsKey("step_id")) throw new IllegalStateException("undefined params.step_id with include+addParams");
if(!params.containsKey("step_name")) throw new IllegalStateException("undefined params.step_name with include+addParams");

include {moduleLoad} from '../../modules/utils/functions.nf'
include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'


workflow PIHAT01 {
	take:
		genomeId
		vcf
		samples
		blacklisted_bed
		remove_samples
	main:
		version_ch = Channel.empty();
		to_zip = Channel.empty();
		vcf2contig_ch = VCF_TO_BED([:],vcf)
		version_ch = version_ch.mix(vcf2contig_ch.version)

		ctgvcf_ch = vcf2contig_ch.bed.splitCsv(header: false,sep:'\t',strip:true).
                               map{T->[T[0],file(T[3])]}.
				combine(samples)


		perCtg = PLINK_PER_CONTIG(genomeId, blacklisted_bed, ctgvcf_ch)
		version_ch = version_ch.mix(perCtg.version.collect())


		pihat_ch = MERGE_PIHAT_VCF(vcf, remove_samples, perCtg.vcf.map{T->T.join("\t")}.collect())
		version_ch = version_ch.mix(pihat_ch.version)


		to_multiqc = pihat_ch.images.mix(pihat_ch.multiqc)

		version_ch = MERGE_VERSION("pihat",version_ch.collect())

	emit:
		version = version_ch
		plink_genome = pihat_ch.plink_genome
		to_multiqc = to_multiqc
		genome_bcf = pihat_ch.genome_bcf /* optional, used for PCA */
	}




process PLINK_PER_CONTIG {
tag "${vcf.name}/${contig}"
afterScript "rm -rf TMP"
memory "3g"
input:
	val(genomeId)
	path("BLACKLISTED/*")
        tuple val(contig),path(vcf),path(samples)
output:
	tuple val(contig),path("${contig}.bcf"),emit:vcf
	path("version.xml"),emit:version
when:
	contig.matches(params.pihat.contig_regex)
script:
	def genome = params.genomes[genomeId]
	if(genome==null) throw new IllegalArgumentException("cannot get genome[${genomeId}]")
	def reference = genome.fasta
	def filters = params.pihat.filters
	def pihatmaf = (params.pihat.MAF as double)
	def pihatMinGQ = (params.pihat.min_GQ as int)
	def f_missing= (params.pihat.f_missing as double)
	def minDP= (params.pihat.min_DP as int)
	def maxDP= (params.pihat.max_DP as int)
	def gnomad_genome_path = genome.gnomad_genome
	def gnomad_exome_path =  genome.gnomad_exome

if(contig.matches("(chr)?[0-9]+"))

"""
hostname 1>&2
${moduleLoad("plink bcftools jvarkit bedtools")}

set -x

mkdir TMP

function countIt {
	echo -n "COUNT : " 1>&2
	bcftools query -N -f a "\$1" | wc -c 1>&2
	}


# reduce size of exclude bed
if test ! -L "BLACKLISTED/NO_FILE" ; then

	awk -F '\t' '(\$1=="${contig}")' BLACKLISTED/*.bed > TMP/jeter.x.bed

else

	touch TMP/jeter.x.bed

fi

# create fake file if empty
if [ ! -s TMP/jeter.x.bed ] ; then

	echo "${contig}\t0\t1" > TMP/jeter.x.bed

fi

bcftools view -m2 -M2 ${filters} --types snps -O u \
		${samples.name.equals("NO_FILE")?"":"--samples-file '${samples.toRealPath()}'"} \
		"${vcf.toRealPath()}" "${contig}" |\
	bcftools view --targets-file ^TMP/jeter.x.bed  --targets-overlap 2 --exclude-uncalled  --min-af "${pihatmaf}" --max-af "${1.0 - (pihatmaf as Double)}"  -i 'AC>0 ${contig.matches("(chr)?Y")?"":"&& F_MISSING < ${f_missing}"}'  -O v |\
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/vcffilterjdk.jar --nocode -e 'final double dp= variant.getGenotypes().stream().filter(G->G.isCalled() && G.hasDP()).mapToInt(G->G.getDP()).average().orElse(${minDP}); if( dp<${minDP} || dp>${maxDP}) return false; if (variant.getGenotypes().stream().filter(G->G.isCalled() && G.hasGQ()).anyMatch(G->G.getGQ()< ${pihatMinGQ} )) return false; return true;' |\
	bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -x '^INFO/AC,INFO/AF,INFO/AN,QUAL,^FORMAT/GT' -O b -o TMP/jeter2.bcf
	mv -v TMP/jeter2.bcf TMP/jeter1.bcf

	countIt TMP/jeter1.bcf

	for G in "${gnomad_genome_path}" "${gnomad_exome_path}"
	do
	if [ -f "\${G}"  ] ; then

		echo "##\${G}"

		# extract bed for this vcf, extends to avoid too many regions
		bcftools query -f '%CHROM\t%POS0\t%END\\n' TMP/jeter1.bcf |\
			java -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "\${G}" --column 1 --convert SKIP |\
			sort -T TMP -t '\t' -k1,1 -k2,2n |\
			bedtools merge -i - -d 1000 > TMP/gnomad.bed
	
		if [ ! -s gnomad.bed ] ; then
        		echo "${contig}\t0\t1" > TMP/gnomad.bed
		fi

		# get filtered variants in gnomad
		bcftools query -e 'FILTER=="." || FILTER=="PASS"'  --regions-file TMP/gnomad.bed -f '%CHROM\t%POS0\t%END\\n' "\${G}" |\
			java  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP |\
			sort -T	TMP -t '\t' -k1,1	-k2,2n > TMP/x.gnomad.bed
		
		if [ ! -s x.gnomad.bed ] ; then
        		echo "${contig}\t0\t1" > TMP/x.gnomad.bed
		fi
		
		bcftools view --targets-file "^TMP/x.gnomad.bed" --targets-overlap 2 -O b -o TMP/jeter2.bcf TMP/jeter1.bcf
		mv TMP/jeter2.bcf TMP/jeter1.bcf
	
		countIt TMP/jeter1.bcf
	fi
	done
	
	bcftools view -O b -o TMP/jeter.bcf TMP/jeter1.bcf
	rm TMP/jeter1.bcf

	countIt TMP/jeter.bcf

#
# Dans le fichier hardyweinberg.hwe enlever les variants dont la colonne P < 0.00001
#
plink --double-id --bcf TMP/jeter.bcf --allow-extra-chr --hardy --allow-no-sex --out "TMP/hardyweinberg.txt"

awk '(\$9 <  0.00001) {print \$2}' "TMP/hardyweinberg.txt.hwe" > TMP/xclude_ids.txt

if [ -s "TMP/xclude_ids.txt" ] ; then
	bcftools view -e 'ID=@TMP/xclude_ids.txt' -O b -o TMP/jeter2.bcf TMP/jeter.bcf
	mv TMP/jeter2.bcf TMP/jeter.bcf
	countIt TMP/jeter.bcf
fi

#
# sélection de SNP indépendants 
#
plink --double-id --bcf  TMP/jeter.bcf --allow-extra-chr  --indep-pairwise 50 10 0.2 --out TMP/plink


plink --double-id --bcf  TMP/jeter.bcf --allow-extra-chr  --extract TMP/plink.prune.in --recode vcf --out TMP/pruned_data
test -f TMP/pruned_data.vcf

plink --double-id --bcf  TMP/jeter.bcf --allow-extra-chr --r2 --ld-window 50 --ld-window-kb 5000 --ld-window-r2 0.2 --out TMP/ld
test -f TMP/ld.ld

awk '{print \$3;}'  TMP/ld.ld  | uniq  | sort -T . | uniq > TMP/snpInLD.txt

plink --double-id --vcf TMP/pruned_data.vcf --exclude TMP/snpInLD.txt --recode vcf --out TMP/indepSNP_data

bcftools view -O b -o "${contig}.bcf" TMP/indepSNP_data.vcf
bcftools index "${contig}.bcf"
countIt "${contig}.bcf"


cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">prepare pihat data per autosome</entry>
	<entry key="vcf">${vcf}</entry>
	<entry key="samples">${samples.toRealPath()}</entry>
	<entry key="contig">${contig}</entry>
	<entry key="maf">${pihatmaf}</entry>
	<entry key="min.GQ">${pihatMinGQ}</entry>
	<entry key="F_MISSING">${f_missing}</entry>
	<entry key="min.DP">${minDP}</entry>
	<entry key="max.DP">${maxDP}</entry>
	<entry key="filters">${filters}</entry>
</properties>
EOF

"""

else if(contig.matches("(chr)?[XY]"))

"""
hostname 1>&2
${moduleLoad("bcftools")}

mkdir TMP

# reduce size of exclude bed

if ${!blacklisted.name.equals("NO_FILE")}  ; then
awk -F '\t' '(\$1=="${contig}")' "${blacklisted}" > TMP/jeter.x.bed
fi

if [ ! -s TMP/jeter.x.bed ] ; then
        echo "${contig}\t0\t1" > TMP/jeter.x.bed
fi


bcftools view -m2 -M2 --apply-filters '.,PASS' --types snps -O u \
		${samples.name.equals("NO_FILE")?"":"--samples-file '${samples.toRealPath()}'"} \
		"${vcf.toRealPath()}" "${contig}" |\
	bcftools view --targets-file ^TMP/jeter.x.bed  --targets-overlap 1 --exclude-uncalled  -i 'AC>0' -O u |\
	bcftools annotate -x 'INFO,ID,FILTER,QUAL,^FORMAT/GT' -O b -o TMP/jeter.bcf
	

mv -v TMP/jeter.bcf "${contig}.bcf"
bcftools index "${contig}.bcf"


cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">prepare pihat data per sexual contig</entry>
	<entry key="vcf">${vcf}</entry>
	<entry key="samples">${samples}</entry>
	<entry key="contig">${contig}</entry>
</properties>
EOF
"""
else
"""
touch "${contig}.bcf"
touch "${contig}.bcf.csi"

echo "<properties/>" > version.xml
"""
stub:
"""
touch "${contig}.bcf" "${contig}.bcf.csi" version.xml
"""
}



process MERGE_PIHAT_VCF {
tag "N=${L.size()}"
cache "lenient"
afterScript ""
input:
	path(vcf)
	path(remove_samples)
	val(L)
output:
	path("${prefix}plink.genome.txt.gz"),emit:plink_genome
	path("genome.bcf"),optional:true,emit:genome_bcf
	//multiqc 
	path("${prefix}*.png"),emit:images
        path("${prefix}multiqc.*"),emit:multiqc
	path("version.xml"),emit:version
	        
script:
	prefix = (params.prefix?:"")+params.step_id+"."

	maxPiHat = (params.pihat.pihat_max as double)
	def save_genome_vcf = task.ext.save_genome_vcf

	def whatispihat = "The probable relatives and duplicates are detected based on pairwise identify-by-state (IBS) from which a variable called PIHAT is calculated via PLINK"

if(!L.isEmpty())
"""
hostname 2>&1
${moduleLoad("plink bcftools r/3.6.3")}
set -x

mkdir -p TMP

cat << EOF > TMP/jeter.list
${L.join("\n")}
EOF

awk -F '\t' '(\$1 ~ /^(chr)?[0-9]+\$/ ) {print \$2}' TMP/jeter.list > TMP/autosomes.list 

# concat autosomes
bcftools concat --allow-overlaps  -O b --file-list TMP/autosomes.list -o "TMP/jeter.bcf"
plink --double-id --bcf TMP/jeter.bcf  --allow-extra-chr --genome --out TMP/plink ${remove_samples.name.equals("NO_FILE")?"":"--remove \"${remove_samples}\""} 


## for ACP subworkflow, I may need the vcf that was used to create the genome
if ${save_genome_vcf} ; then
	cp -v TMP/jeter.bcf genome.bcf
	bcftools index -f genome.bcf
fi


# extract VCF for sexual chromsomes
awk -F '\t' '(\$1 ~ /^(chr)?[XY]+\$/ ) {print \$2}' TMP/jeter.list > TMP/sex.list

if [ -s "TMP/sex.list" ] ; then
	bcftools concat --allow-overlaps  -O b --file-list TMP/sex.list -o "TMP/jeter.bcf"
	plink --double-id --bcf TMP/jeter.bcf  --allow-extra-chr --make-bed --impute-sex --out TMP/impute-sex
fi

## plot it
awk '(NR>1) {printf("%s_%s\\t%s\\n",\$1,\$3,\$10);}' TMP/plink.genome > TMP/jeter.tsv



# create table sample/avg(pihat)/status
awk '(NR>1) {P[\$1]+=1.0*(\$10);P[\$3]+=1.0*(\$10);C[\$1]++;C[\$3]++;} END{for(S in P) printf("%s\t%f\\n",S,P[S]/C[S]);}' TMP/plink.genome |\\
	LC_ALL=C sort -T . -t '\t' -k2,2gr > "TMP/sample2avg.pihat.tsv"

# create table with KEEP / DISCARD status
awk -F '\t' 'function sn(SN) {nuscore=split(SN,a,/_/); nuscore=nuscore/4; sample="";for(i=1;i<=nuscore;i++) sample=sprintf("%s%s%s",sample,(i==1?"":"_"),a[i]); return sample;}  {printf("%s\t%s\t%s\\n",sn(\$1),\$2,(\$2 < ${maxPiHat} ? "KEEP":"DISCARD"));}' TMP/sample2avg.pihat.tsv > TMP/samples.keep.status

# the following test will fail if all samples are above maxPiHat
cut -f 3 TMP/samples.keep.status | grep -F -w KEEP -m1 1>&2

cat << EOF > TMP/jeter.R
png("${prefix}plot.pihat.png")
genome <- read.table(file="TMP/jeter.tsv",sep="\\t",header=FALSE)
plot(genome\\\$V2,ylim=c(0,1.0),xlab="Individuals Pair", ylab="PI-HAT", main="${prefix}PI-HAT")
abline(h=${maxPiHat},col="blue");
dev.off()

T1<-read.table("TMP/sample2avg.pihat.tsv",sep="\\t",header=FALSE,col.names=c("S","X"),colClasses=c("character","numeric"))
head(T1)
png("${prefix}plot.sample2avg.pihat.png")
boxplot(T1\\\$X ,ylim=c(0,max(T1\\\$X)),main="${prefix}AVG(PIHAT)/SAMPLE",sub="${vcf.name}",xlab="Sample",ylab="pihat")
abline(h=${maxPiHat},col="blue");
dev.off()

EOF

R --no-save < TMP/jeter.R




##
## create MULTIQC CONFIG
##
cat << EOF > "${prefix}multiqc.config.yaml"
custom_data:
  pihat_plot1_${params.step_id}:
    parent_id: pihat_section_${params.step_id}
    parent_name: "PIHAT ${params.step_name}"
    parent_description: "${whatispihat}"
    section_name: "Pihat  ${params.step_name}"
    description: "plot of pihat"
  pihat_plot2_${params.step_id}:
    parent_id: pihat_section_${params.step_id}
    parent_name: "PIHAT  ${params.step_name}"
    parent_description: "${whatispihat}"
    section_name: "Sample to Average pihat"
    description: "Sample to Average pihat"
sp:
  pihat_plot1_${params.step_id}:
    fn: "${prefix}plot.pihat.png"
  pihat_plot2_${params.step_id}:
    fn: "${prefix}plot.sample2avg.pihat.png"
ignore_images: false
EOF



cat << EOF > TMP/${prefix}multiqc.genome_mqc.html
<!--
id: 'highpihat_${params.step_id}'
parent_id: pihat_section_${params.step_id}
parent_name: "PIHAT ${params.step_name}"
parent_description: "${whatispihat}"
section_name: 'High Pihat ${params.step_name}'
description: ' 10 higher pihat pairs.'
-->
EOF

awk '(NR>1) {printf("%s\t%s\t%s\\n",\$1,\$3,\$10);}' TMP/plink.genome |\
	LC_ALL=C sort -t '\t' -T TMP -k3,3gr |\
	head -n 10 |\
	awk -F '\t' 'function sn(SN) {nuscore=split(SN,a,/_/); nuscore = nuscore/4; sample="";for(i=1;i<=nuscore;i++) sample=sprintf("%s%s%s",sample,(i==1?"":"_"),a[i]); return sample;} BEGIN {printf("<table class=\\"table\\"><tr><th>SN1</th><th>SN2</th><th>PIHAT</th></tr>\\n");} {printf("<tr><td>%s</td><td>%s</td><td>%s</td></tr>\\n",sn(\$1),sn(\$2),\$3);} END {printf("</table>\\n");}' >> TMP/${prefix}multiqc.genome_mqc.html


cat << EOF > TMP/${prefix}multiqc.removed_samples_mqc.html
<!--
id: 'rmsamples_${params.step_id}'
parent_id: pihat_section_${params.step_id}
parent_name: "PIHAT ${params.step_name}"
parent_description: "${whatispihat}"
section_name: 'Removed Samples ${params.step_name}'
description: 'Samples that would be removed with max.pihat=${params.pihat.pihat_max}.'
-->
EOF

awk -F '\t' 'BEGIN {printf("<table class=\\"table\\"><tr><th>Sample</th><th>Pihat</th></tr>\\n");} (\$3=="DISCARD") {printf("<tr><td>%s</td><td>%s</td></tr>\\n",\$1,\$2);} END {printf("</table>\\n");}' TMP/samples.keep.status  >> TMP/${prefix}multiqc.removed_samples_mqc.html


gzip --best TMP/plink.genome
mv TMP/plink.genome.gz "${prefix}plink.genome.txt.gz"

mv TMP/${prefix}multiqc.genome_mqc.html ./
mv TMP/${prefix}multiqc.removed_samples_mqc.html ./

###########################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry id="name">${task.process}</entry>
	<entry key="description">merge pihat data per contig, create pihat data</entry>
</properties>
EOF
"""

stub:
"""
touch "${prefix}plink.genome.txt.gz" "${prefix}pihat.png" "${prefix}removed_samples.txt"
echo "<properties/>" > version.xml
"""
}

