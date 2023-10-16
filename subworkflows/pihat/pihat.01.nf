/*

Copyright (c) 2023 Pierre Lindenbaum

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


include {moduleLoad} from '../../modules/utils/functions.nf'
include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'


workflow PIHAT01 {
	take:
		genomeId
		vcf
		samples
	main:
		version_ch = Channel.empty();
		to_zip = Channel.empty();
		vcf2contig_ch = VCF_TO_BED([with_header:false],vcf)
		version_ch = version_ch.mix(vcf2contig_ch.version)

		ctgvcf_ch = vcf2contig_ch.bed.splitCsv(header: false,sep:'\t',strip:true).
                               map{T->[T[0],file(T[3])]}.
				combine(samples)


		perCtg = PLINK_PER_CONTIG(genomeId, ctgvcf_ch)
		version_ch = version_ch.mix(perCtg.version.collect())


		pihat_ch = MERGE_PIHAT_VCF(vcf, perCtg.vcf.map{T->T.join("\t")}.collect())
		version_ch = version_ch.mix(pihat_ch.version)

		mqc_ch = MULTIQC([:],vcf,pihat_ch.pihat_png, pihat_ch.pihat_sample2avg_png, pihat_ch.pihat_removed_samples, pihat_ch.plink_genome)
		version_ch = version_ch.mix(mqc_ch.version)

		version_ch = MERGE_VERSION("pihat",version_ch.collect())

	emit:
		version = version_ch
		pihat_png = pihat_ch.pihat_png
		pihat_removed_samples = pihat_ch.pihat_removed_samples
		plink_genome = pihat_ch.plink_genome
		pihat_sample2avg_png  = pihat_ch.pihat_sample2avg_png
		pihat_multiqc_zip = mqc_ch.zip
	}




process PLINK_PER_CONTIG {
tag "${vcf.name}/${contig}"
afterScript "rm -rf TMP"
memory "3g"
input:
	val(genomeId)
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
	def blacklisted = params.pihat.blacklisted
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
	bcftools query -f a "\$1" | wc -c 1>&2
	}


# reduce size of exclude bed
if [ -f "${blacklisted}" ] ; then

	awk -F '\t' '(\$1=="${contig}")' "${blacklisted}" > TMP/jeter.x.bed

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

if [ -s "${blacklisted}" ] ; then
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
	val(L)
output:
	path("${prefix}pihat.png"),emit:pihat_png
	path("${prefix}sample2avg.pihat.png"),emit:pihat_sample2avg_png
	path("${prefix}removed_samples.txt"),emit:pihat_removed_samples
	path("${prefix}plink.genome.txt.gz"),emit:plink_genome
	path("version.xml"),emit:version
script:
	prefix = params.prefix?:""
	maxPiHat = (params.pihat.pihat_max as double)

if(!L.isEmpty())
"""
hostname 2>&1
${moduleLoad("plink bcftools r/3.6.3")}
set -x

mkdir TMP

cat << EOF > TMP/jeter.list
${L.join("\n")}
EOF

awk -F '\t' '(\$1 ~ /^(chr)?[0-9]+\$/ ) {print \$2}' TMP/jeter.list > TMP/autosomes.list 

bcftools concat --allow-overlaps  -O b --file-list TMP/autosomes.list -o "TMP/jeter.bcf"
plink --double-id --bcf TMP/jeter.bcf  --allow-extra-chr --genome --out TMP/plink

awk -F '\t' '(\$1 ~ /^(chr)?[XY]+\$/ ) {print \$2}' TMP/jeter.list > TMP/sex.list

if [ -s "TMP/sex.list" ] ; then
	bcftools concat --allow-overlaps  -O b --file-list TMP/sex.list -o "TMP/jeter.bcf"
	plink --double-id --bcf TMP/jeter.bcf  --allow-extra-chr --make-bed --impute-sex --out TMP/impute-sex
fi

## plot it
awk '(NR>1) {printf("%s_%s\\t%s\\n",\$1,\$3,\$10);}' TMP/plink.genome > TMP/jeter.tsv



# create table sample/avg(pihat)/status
awk '(NR>1) {P[\$1]+=1.0*(\$10);P[\$3]+=1.0*(\$10);C[\$1]++;C[\$3]++;} END{for(S in P) printf("%s\t%f\\n",S,P[S]/C[S]);}' TMP/plink.genome |\
	LC_ALL=C sort -T . -t '\t' -k2,2gr > "${prefix}sample2avg.pihat.tsv"


awk -F '\t' '(\$2 < ${maxPiHat}) {print \$1;}' "${prefix}sample2avg.pihat.tsv" | sort -T TMP > TMP/jeter.keep.samples.txt
awk -F '\t' '(\$2 >= ${maxPiHat})' "${prefix}sample2avg.pihat.tsv" > "${prefix}removed_samples.txt"

# the following test will fail if all samples are above maxPiHat
test -s TMP/jeter.keep.samples.txt

cat << EOF > TMP/jeter.R
png("${prefix}pihat.png")
genome <- read.table(file="TMP/jeter.tsv",sep="\\t",header=FALSE)
plot(genome\\\$V2,ylim=c(0,1.0),xlab="Individuals Pair", ylab="PI-HAT", main="${prefix}PI-HAT")
abline(h=${maxPiHat},col="blue");
dev.off()

T1<-read.table("${prefix}sample2avg.pihat.tsv",sep="\\t",header=FALSE,col.names=c("S","X"),colClasses=c("character","numeric"))
head(T1)
png("${prefix}sample2avg.pihat.png")
boxplot(T1\\\$X ,ylim=c(0,max(T1\\\$X)),main="${prefix}AVG(PIHAT)/SAMPLE",sub="${vcf.name}",xlab="Sample",ylab="pihat")
abline(h=${maxPiHat},col="blue");
dev.off()

EOF

R --no-save < TMP/jeter.R


gzip --best TMP/plink.genome
mv TMP/plink.genome.gz "${prefix}plink.genome.txt.gz"



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

process MULTIQC {
input:
	val(meta)
	path(vcf)
        path(pihat_png)
        path(pihat_sample2avg_png)
        path(pihat_removed_samples)
	path(genome_tsv_gz)
output:
	path("${params.prefix?:""}multiqc.zip"),emit:zip
	path("version.xml"),emit:version
script:
	def m = 20
	def prefix =params.prefix?:""
	def prefix0 = prefix.endsWith(".")?prefix.substring(0,prefix.length()-1):prefix
"""
${moduleLoad("multiqc")}
hostname 1>&2
mkdir -p TMP

cat << EOF > multiqc_config.yaml
custom_data:
  pihat01:
    section_name: "Pihat"
    description: "Pihat"
  pihat02:
    section_name: "Sample to Average Pihat"
    description: "Sample to Average Pihat"
sp:
  pihat01:
    fn: "${pihat_png.name}"
  pihat02:
    fn: "${pihat_sample2avg_png.name}"
ignore_images: false
EOF


cat << EOF > TMP/genome_mqc.html
<!--
id: 'mytable'
section_name: 'High Pihat'
description: '${m} higher pihat pairs.'
-->
EOF



gunzip -c '${genome_tsv_gz}' |\
	awk '(NR>1) {printf("%s\t%s\t%s\\n",\$1,\$3,\$10);}' |\
	LC_ALL=C sort -t '\t' -T TMP -k3,3gr |\
	head -n ${m} |\
	awk -F '\t' 'BEGIN {printf("<table><tr><th>SN1</th><th>SN2</th><th>PIHAT</th></tr>\\n");} {printf("<tr><td>%s</td><td>%s</td><td>%s</td></tr>\\n",\$1,\$2,\$3);} END {printf("</table>\\n");}' >> TMP/genome_mqc.html


cat << EOF > TMP/removed_samples_mqc.html
<!--
id: 'rmsamples'
section_name: 'Removed Samples'
description: 'Samples that would be removed with max.pihat=${params.pihat.pihat_max}.'
-->
EOF

awk -F '\t' 'BEGIN {printf("<ol>\\n");} {printf("<li>%s</li>\\n",\$0);} END {printf("</ol>\\n");}' '${pihat_removed_samples}'  >> TMP/removed_samples_mqc.html

mkdir -p "${params.prefix}multiqc"

multiqc   --outdir "${params.prefix}multiqc"  --filename "${prefix0}.multiqc"  --title "${params.prefix}PIHAT for ${vcf.name}" \
	--comment "  PLINK --genome estimates relatedness of all pairs of samples and reports identify by decent (IBD, a measure of whether identical regions of two genomes were inherited from the same ancestry) in the PI_HAT (actually, proportional IBD, i.e. P(IBD=2) + 0.5*P(IBD=1)) column of the result file. A PI_HAT value close to 1 would indicate a duplicate sample. VCF was ${vcf}." \
	--interactive --force --verbose --module custom_content .

zip -9 -r "${params.prefix}multiqc.zip" "${params.prefix}multiqc"

##################################################################################
cat <<- EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">MultiQC</entry>
        <entry key="wget.version">\$( multiqc --version )</entry>
</properties>
EOF
"""
}
