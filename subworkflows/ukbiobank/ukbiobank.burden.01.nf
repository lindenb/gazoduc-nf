include {moduleLoad;parseBoolean} from '../../modules/utils/functions.nf'
include {UKBIOBANK_VCF2BED_01} from '../../modules/ukbiobank/ukbiobank.vcf2bed.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {RVTESTS_POST_PROCESS} from '../rvtests/rvtests.post.process.01.nf'
include {CONCAT_FILES_01 as CONCAT1; CONCAT_FILES_01 as CONCAT2} from '../../modules/utils/concat.files.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'




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

		each_gene = genes_ch.output.splitCsv(header:false,sep:'\t')

		ukvcf2bed_ch = UKBIOBANK_VCF2BED_01(meta)
		version_ch = version_ch.mix(ukvcf2bed_ch.version)

		test_gene_ch = TEST_ONE_GENE(meta, reference, gff3 ,pedigree, ukvcf2bed_ch.bed, each_gene)
		version_ch = version_ch.mix(test_gene_ch.version)

		concat_ch = CONCAT1(meta,test_gene_ch.output.collect())
		version_ch = version_ch.mix(concat_ch.version)

		variants_ch = CONCAT2(meta,test_gene_ch.variants.collect())
		version_ch = version_ch.mix(variants_ch.version)

		digest_ch = RVTESTS_POST_PROCESS(meta, reference,"UKVCF" ,concat_ch.output)
                version_ch = version_ch.mix(digest_ch.version)

		version_ch = MERGE_VERSION(meta, "burden UK", "Burden UK", version_ch.collect())
		
		html = VERSION_TO_HTML(params,version_ch.version)


		to_zip = Channel.empty().mix(version_ch).
			mix(variants_ch.output).
			mix(digest_ch.zip)
			mix(html.html)

	emit:
		version = version_ch
		zip = to_zip
	}


process EXTRACT_GENES {
	input:
		val(meta)
		val(reference)
		path(gff)
	output:
		path("genes.tsv"),emit:output
		path("version.xml"),emit:version
	script:
	"""
	hostname 1>&2
	${moduleLoad("jvarkit bedtools")}
	set -o pipefail

	java -jar \${JVARKIT_DIST}/gtf2bed.jar "${gff}" -c 'biotype,Parent,Name,transcript_id' |\
		awk -F '\t' '(\$4=="protein_coding" && \$5 ~ /^gene\\:/ && \$7!=".") {OFS="\t";if(\$6=="." || \$6=="") \$6=\$5;print;}' |\
		sed 's/\tgene:/\t/g' |\
		sort -T . -t '\t' -k1,1 -k5,5 |\
		bedtools groupby -g 1,5 -c 2,3,6,7  -o min,max,first,distinct > genes.tsv

	test -s genes.tsv
	cat <<- EOF > version.xml
	<properties id="${task.process}">
        	<entry key="name">${task.process}</entry>
        	<entry key="description">extract genes and transcripts</entry>
		<entry key="gff3">${gff}</entry>
		<entry key="count">\$(wc -l < genes.tsv)</entry>
	</properties>
	EOF
	"""
	}

process TEST_ONE_GENE {
tag "${contig}:${start}-${end} ${ensg}/${geneName}"
afterScript "rm -rf TMP"
errorStrategy "retry"
maxRetries 5
memory "5g"
input:
	val(meta)
	val(reference)
	path(gff3)
	path(pedigree)
	path(ukbiobank2bed)
	tuple val(contig),val(ensg),val(start),val(end),val(geneName),val(enst)

output:
	path("assoc.list"),emit:output
	path("variants.bed"),emit:variants
	path("version.xml"),emit:version
when:
	!ensg.equals("ENSG00000155657")
script:
	def soacn=meta.soacn?:"SO:0001574,SO:0001575,SO:0001818"
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
	# SLS: a negative person ID in the FAM file means that the corresponding participant has withdrawn consent and should therefore be excluded.
	bcftools query -l "\${F}" | awk '!(\$1 ~ /^\-/)' | sort -T TMP | uniq > TMP/a
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

# case/list for contrast
awk -F '\t' '(\$6=="1") {printf("%s\\n",\$2);}' TMP/jeter.ped > TMP/ctrl.list
awk -F '\t' '(\$6=="2") {printf("%s\\n",\$2);}' TMP/jeter.ped > TMP/cases.list

# problem with uk biobank is we cannot use bcftools concat because they don't have the same samples !!
i=1
cut -f 4 TMP/vcf.bed | while read F
do

	bcftools view -O u --regions-file TMP/jeter.bed --samples-file TMP/samples.txt --force-samples "\${F}" |\
		bcftools view -i 'AC>0 && F_MISSING<0.01' --trim-alt-alleles --max-af '0.01' -O b -o TMP/chunk.\${i}.bcf
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
bcftools query -N -f '.' TMP/jeter.bcf | wc -c 1>&2




## extract sites only, annotate
echo "${enst}" | tr "," "\\n" | sort | uniq > TMP/each_tr.txt
test -s TMP/each_tr.txt


bcftools csq -O u --force --local-csq --ncsq 10000 --fasta-ref "${reference}" --gff-annot "${gff3}" TMP/jeter.bcf |\
	    bcftools +contrast -0 TMP/ctrl.list -1 TMP/cases.list -a PASSOC,FASSOC,NASSOC,NOVELAL,NOVELGT -O z -o TMP/jeter2.vcf.gz -

bcftools index -t TMP/jeter2.vcf.gz

echo -n "NUMBER OF VARIANTS AFTER ANNOT:" 1>&2
bcftools query -N -f '.' TMP/jeter2.vcf.gz |wc -c 1>&2

# loop over each TRanscript
i=1
cat TMP/each_tr.txt | while read TR
do
	bcftools view TMP/jeter2.vcf.gz |\
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/vcffilterso.jar \
			--remove-attribute  --rmnoatt \
			--acn "${soacn}" |\
	java -jar ${JVARKIT_DIST}/vcfburdenfiltergenes.jar -a "\${TR}" |\
	bcftools view -O z -o TMP/jeter3.vcf.gz

	bcftools index -t -f TMP/jeter3.vcf.gz


	echo "${ensg}_${geneName}_\${TR}\t${contig}:${start}-${end}" > TMP/variants.setfile


	rvtest  --noweb \
	        --inVcf TMP/jeter3.vcf.gz \
		--setFile TMP/variants.setfile \
		--pheno TMP/jeter.ped \
		--out "ASSOC/part.\${i}" \
	        --burden cmc,zeggini,mb,fp,exactCMC,cmcWald,rarecover,cmat \
		--vt price,analytic \
		--kernel 'skat[nPerm=1000],kbac,skato' 1>&2 2> TMP/last.rvtest.log

	echo -n "${contig}\t${start}\t${end}\t${ensg}\t${geneName}\t\${TR}\t" >> variants.bed
	bcftools query -f '%POS:%REF:%ALT:%INFO/PASSOC:%INFO/FASSOC\\n' TMP/jeter3.vcf.gz | paste -s -d ',' >> variants.bed

	i=\$((i+1))
done

find \${PWD}/ASSOC -type f -name "part*assoc" > assoc.list


###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">invoke rvtest for gene</entry>
	<entry key="rvtest.path">\$(which rvtest)</entry>
	<entry key="gene">${ensg}</entry>
	<entry key="geneName">${geneName}</entry>
	<entry key="transcripts">${enst}</entry>
	<entry key="soacn">${soacn}</entry>
	<entry key="rvtest.version">\$(rvtest  --version 2> /dev/null  | head -n1)</entry>
	<entry key="bedtools.version">\$(  bedtools --version )</entry>
	<entry key="uk.bed">${ukbiobank2bed}</entry>
</properties>
EOF
"""
}
