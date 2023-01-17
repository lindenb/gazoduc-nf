/*

Copyright (c) 2022 Pierre Lindenbaum

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


include {moduleLoad;getGnomadExomePath;getGnomadGenomePath} from '../../modules/utils/functions.nf'
include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'

def gazoduc = gazoduc.Gazoduc.getInstance(params).putGnomad()


workflow PIHAT01 {
	take:
		meta
		reference
		vcf
		samples
	main:
		version_ch = Channel.empty();
		to_zip = Channel.empty();
		vcf2contig_ch = VCF_TO_BED(meta,vcf)
		version_ch = version_ch.mix(vcf2contig_ch.version)

		ctgvcf_ch = vcf2contig_ch.bed.splitCsv(header: false,sep:'\t',strip:true).
                               map{T->[T[0],file(T[3])]}.combine(samples)


		perCtg = PLINK_PER_CONTIG(meta,reference, ctgvcf_ch)
		version_ch = version_ch.mix(perCtg.version.collect())


		pihat_ch = MERGE_PIHAT_VCF(meta.plus(["title":"${vcf.name}"]), perCtg.vcf.map{T->T.join("\t")}.collect())
		version_ch = version_ch.mix(pihat_ch.version)

		version_ch = MERGE_VERSION(meta, "pihat", "pihat on ${vcf}", version_ch.collect())

	emit:
		version = version_ch
		pihat_png = pihat_ch.pihat_png
		pihat_removed_samples = pihat_ch.pihat_removed_samples
		plink_genome = pihat_ch.plink_genome
		pihat_pdf  = pihat_ch.pihat_pdf
	}


gazoduc.build("pihat_filters", " --apply-filters '.,PASS' ").
	description("pihat filter").
	menu("pihat").
	put()

gazoduc.build("pihat_MAF", 0.1).
	description("pihat MAF").
	menu("pihat").
    setDouble().
	put()

gazoduc.build("pihat_min_GQ", 20).
	description("pihat min GQ").
	menu("pihat").
    setInt().
	put()

gazoduc.build("pihat_f_missing", 0.05).
	description("pihat fraction of missing genotypes").
	menu("pihat").
    setDouble().
	put()

gazoduc.build("pihat_min_DP", 10).
	description("min Depth").
	menu("pihat").
    setInt().
	put()

gazoduc.build("pihat_max_DP", 300).
	description("max Depth").
	menu("pihat").
    setInt().
	put()


gazoduc.build("pihat_max", 0.05).
	description("max pihat ").
	menu("pihat").
    setDouble().
	put()


process PLINK_PER_CONTIG {
tag "${vcf.name}/${contig}"
afterScript "rm -rf TMP"
memory "3g"
input:
	val(meta)
	val(reference)
        tuple val(contig),path(vcf),path(samples)
output:
	tuple val(contig),path("${contig}.bcf"),emit:vcf
	path("version.xml"),emit:version
when:
	contig.matches(meta.pihat_contig_regex?:"(chr)?[0-9XY]+")
script:

	def filters = meta.pihat_filters
	def blacklisted = meta.blacklisted?:""	
	def pihatmaf = (meta.pihat_MAF as double)
	def pihatMinGQ = (meta.pihat_min_GQ as int)
	def f_missing= (meta.pihat_f_missing as double)
	def minDP= (meta.pihat_min_DP as int)
	def maxDP= (meta.pihat_max_DP as int)
	def gnomad_genome_path = getGnomadExomePath(meta,reference)
	def gnomad_exome_path =   getGnomadGenomePath(meta,reference)

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
plink --bcf TMP/jeter.bcf --allow-extra-chr --hardy --allow-no-sex --out "TMP/hardyweinberg.txt"

awk '(\$9 <  0.00001) {print \$2}' "TMP/hardyweinberg.txt.hwe" > TMP/xclude_ids.txt

if [ -s "TMP/xclude_ids.txt" ] ; then
	bcftools view -e 'ID=@TMP/xclude_ids.txt' -O b -o TMP/jeter2.bcf TMP/jeter.bcf
	mv TMP/jeter2.bcf TMP/jeter.bcf
	countIt TMP/jeter.bcf
fi

#
# sélection de SNP indépendants 
#
plink --bcf  TMP/jeter.bcf --allow-extra-chr  --indep-pairwise 50 10 0.2 --out TMP/plink


plink --bcf  TMP/jeter.bcf --allow-extra-chr  --extract TMP/plink.prune.in --recode vcf --out TMP/pruned_data
test -f TMP/pruned_data.vcf

plink --bcf  TMP/jeter.bcf --allow-extra-chr --r2 --ld-window 50 --ld-window-kb 5000 --ld-window-r2 0.2 --out TMP/ld
test -f TMP/ld.ld

awk '{print \$3;}'  TMP/ld.ld  | uniq  | sort -T . | uniq > TMP/snpInLD.txt

plink --vcf TMP/pruned_data.vcf --exclude TMP/snpInLD.txt --recode vcf --out TMP/indepSNP_data

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
	val(meta)
	val(L)
output:
	path("${prefix}pihat.png"),emit:pihat_png
	path("${prefix}sample2avg.pihat.pdf"),emit:pihat_pdf
	path("${prefix}removed_samples.txt"),emit:pihat_removed_samples
	path("${prefix}plink.genome.txt.gz"),emit:plink_genome
	path("version.xml"),emit:version
script:
	prefix = meta.prefix?:""
	maxPiHat = (meta.pihat_max as double)

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
plink --bcf TMP/jeter.bcf  --allow-extra-chr --genome --out TMP/plink

awk -F '\t' '(\$1 ~ /^(chr)?[XY]+\$/ ) {print \$2}' TMP/jeter.list > TMP/sex.list

if [ -s "TMP/sex.list" ] ; then
	bcftools concat --allow-overlaps  -O b --file-list TMP/sex.list -o "TMP/jeter.bcf"
	plink --bcf TMP/jeter.bcf  --allow-extra-chr --make-bed --impute-sex --out TMP/impute-sex
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
pdf("${prefix}sample2avg.pihat.pdf")
boxplot(T1\\\$X ,ylim=c(0,max(T1\\\$X)),main="${prefix}AVG(PIHAT)/SAMPLE",sub="${meta.title}",xlab="Sample",ylab="pihat")
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
