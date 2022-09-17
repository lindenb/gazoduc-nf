include {moduleLoad;parseBoolean} from '../../modules/utils/functions.nf'
include {UKBIOBANK_VCF2BED_01} from '../../modules/ukbiobank/ukbiobank.vcf2bed.01.nf'


workflow UKBIOBANK_BURDEN_01 {
	take:
		meta
		reference
		pedigree
		gff3
	main:
		version_ch = Channel.empty()
		
		genes_ch = EXTRACT_GENES(meta,reference,gff3)
		version_ch = version_ch.mix(genes_ch.version)

		each_gene = genes_ch.bed.splitCsv(header:false,sep:'\t')

		ukvcf2bed_ch = UKBIOBANK_VCF2BED_01(meta)
		version_ch = version_ch.mix(ukvcf2bed_ch.version)

		test_gene_ch = TEST_ONE_GENE(meta, reference, gff3 ,pedigree, ukvcf2bed_ch.bed, each_gene)
		version_ch = version_ch.mix(test_gene_ch.version)
	emit:
		version = version_ch
	}


process EXTRACT_GENES {
	input:
		val(meta)
		val(reference)
		path(gff)
	output:
		path("genes.bed"),emit:bed
		path("version.xml"),emit:version
	script:
	"""
	hostname 1>&2
	${moduleLoad("jvarkit")}
	set -o pipefail

	java -jar \${JVARKIT_DIST}/gtf2bed.jar "${gff}" -c 'biotype,transcript_id'|\
		awk -F '\t' '(\$4=="protein_coding" && \$5!=".")' |\
		cut -f 1,2,3,5 |\
		sort -t '\t' -T . -k1,1 -k2,2 |\
		uniq > genes.bed
	test -s genes.bed
	cat <<- EOF > version.xml
	EOF
	"""
	}

process TEST_ONE_GENE {
tag "${contig}:${start}-${end} ${enst}"
afterScript "rm -rf TMP"
memory "5g"
input:
	val(meta)
	val(reference)
	path(gff3)
	path(pedigree)
	path(ukbiobank2bed)
	tuple val(contig),val(start),val(end),val(enst)
output:
	path("assoc.list"),emit:output
	path("version.xml"),emit:version
script:
	def soacn=meta.soacn?"":"SO:0001574,SO:0001575,SO:0001818"
"""
hostname 1>&2
${moduleLoad("bcftools jvarkit rvtests bedtools")}

set -x

mkdir TMP
mkdir ASSOC

tail -n+2 "${pedigree}" | cut -f 2 | sort -T TMP | uniq > TMP/samples.txt

echo "${contig}\t${start}\t${end}" > TMP/jeter.bed

bedtools intersect -u -a "${ukbiobank2bed}" -b TMP/jeter.bed > TMP/vcf.bed

if test ! -s TMP/vcf.bed ; then

	tail -n1 "${ukbiobank2bed}" > TMP/vcf.bed

fi

# reduce samples
cut -f 4 TMP/vcf.bed | while read F
do
	bcftools query -l "\${F}" | sort -T TMP | uniq > TMP/a
	comm -12 TMP/a TMP/samples.txt |  sort -T TMP | uniq > TMP/b
	mv TMP/b TMP/samples.txt
	test -s TMP/samples.txt
done

test -s TMP/samples.txt

# reformat pedigree
head -n 1 "${pedigree}" > TMP/jeter.ped
tail -n+2 "${pedigree}" |\
	sort -T TMP -t '\t' -k2,2 |\
	join -1 2 -2 1 -t '\t' -o '1.1,1.2,1.3,1.4,1.5,1.6' - TMP/samples.txt >> TMP/jeter.ped

# problem with uk biobank is we cannot use bcftools concat because they don't have the same samples !!
i=1
cut -f 4 TMP/vcf.bed | while read F
do

	bcftools view -O u --regions-file TMP/jeter.bed --samples-file TMP/samples.txt --force-samples "\${F}" |\
		bcftools view -i 'AC>0' --trim-alt-alleles --max-af '0.01' -O b -o TMP/chunk.\${i}.bcf
	bcftools index TMP/chunk.\${i}.bcf
	echo "TMP/chunk.\${i}.bcf" >> TMP/chunk.list

	i=\$((i+1))
done

## merge bcf for this gene

if [ "\${i}" -eq 2 ] ; then

	mv TMP/chunk.1.bcf TMP/jeter.bcf
	mv TMP/chunk.1.bcf.csi TMP/jeter.bcf.csi

else

	bcftools concat -a --file-list TMP/chunk.list -O b -o TMP/jeter.bcf
	bcftools index TMP/jeter.bcf
	rm -fv TMP/chunk.*.bcf TMP/chunk.*.bcf.csi
fi

echo -n "NUMBER OF VARIANTS BEFORE ANNOT:" 1>&2
bcftools query -f '.' TMP/jeter.bcf | wc -c 1>&2


## extract sites only, annotate
bcftools csq -O v --force --local-csq --ncsq 10000 --fasta-ref "${reference}" --gff-annot "${gff3}" TMP/jeter.bcf |\
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/vcffilterso.jar \
			--remove-attribute  --rmnoatt \
			--acn "${soacn}" |\
	java -jar ${JVARKIT_DIST}/vcfburdenfiltergenes.jar -a "${enst}" |\
	bcftools view -O z -o TMP/jeter2.vcf.gz

bcftools index -t TMP/jeter2.vcf.gz

echo -n "NUMBER OF VARIANTS AFTER ANNOT:" 1>&2
bcftools query -f '.' TMP/jeter2.vcf.gz |wc -c 1>&2

echo "${enst}\t${contig}:${start}-${end}" > TMP/variants.setfile


rvtest  --noweb \
        --inVcf TMP/jeter2.vcf.gz \
	--setFile TMP/variants.setfile \
	--pheno TMP/jeter.ped \
	--out "ASSOC/part" \
        --burden cmc,zeggini,mb,fp,exactCMC,cmcWald,rarecover,cmat \
	--vt price,analytic \
	--kernel 'skat[nPerm=1000],kbac,skato' 1>&2 2> TMP/last.rvtest.log

find \${PWD}/ASSOC -type f -name "part*assoc" > assoc.list

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">invoke rvtest for gene</entry>
	<entry key="rvtest.path">\$(which rvtest)</entry>
	<entry key="transcript">${enst}</entry>
	<entry key="soacn">${soacn}</entry>
	<entry key="rvtest.version">\$(rvtest  --version 2> /dev/null  | head -n1)</entry>
	<entry key="bedtools.version">\$(  bedtools --version )</entry>
	<entry key="uk.bed">${ukbiobank2bed}</entry>
</properties>
EOF
"""
}
