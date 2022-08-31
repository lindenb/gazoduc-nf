include {moduleLoad} from '../../modules/utils/functions.nf'

params.reference=""

workflow {
	rsrc_ch = Channel.from([
		[
		"name":"GRCm38",
		"chain":"http://hgdownload.cse.ucsc.edu/goldenpath/mm10/liftOver/mm10ToHg19.over.chain.gz",
		"src_mapping":"https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCm38_ensembl2UCSC.txt",
		"src_inverse_mapping":true,
		"dest_mapping":"https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh37_gencode2UCSC.txt",
		"dest_inverse_mapping":true,
		"vcf":"http://ftp.ensembl.org/pub/release-102/variation/vcf/mus_musculus/mus_musculus.vcf.gz",
		"fasta":"http://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz",
		]
		])
	ch1 = DOWNLOAD_AND_LIFT(params,rsrc_ch.map{T->T.plus("reference":params.reference)})
	ch2 = ch1.output.splitCsv(header:false,sep:'\t')

	ch3 = APPLY_LIFTOVER(params,ch2)
	
	CONCATENATE(params,ch3.output.groupTuple())
	}

process DOWNLOAD_AND_LIFT {
tag "${row.name}"
afterScript "rm -rf TMP"
memory "5g"
input:
	val(meta)
	val(row)
output:
	path("${row.name}.bed"),emit:output
script:
"""
hostname 1>&2
${moduleLoad("samtools jvarkit bcftools picard")}
mkdir -p TMP


wget -q -O TMP/ref.fa.gz "${row.fasta}"

samtools dict -o TMP/ref.dict TMP/ref.fa.gz

rm TMP/ref.fa.gz

wget -O - "${row.src_mapping}" |\
	awk -F '\t' '{if(\$1=="" || \$2=="") next; printf("%s\t%s\\n",${row.src_inverse_mapping?"\$2,\$1":"\$1,\$2"});}' |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f TMP/ref.dict --column 2 --convert SKIP 	> TMP/src.mapping

wget -O - "${row.dest_mapping}" |\
	awk -F '\t' '{if(\$1=="" || \$2=="") next; printf("%s\t%s\\n",${row.dest_inverse_mapping?"\$2,\$1":"\$1,\$2"});}' |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${row.reference}" --column 2 --convert SKIP 	> TMP/dest.mapping


wget -q -O - "${row.chain}" |\
	gunzip -c |\
	java -jar /${JVARKIT_DIST}/convertliftoverchain.jar \
	-R1 TMP/src.mapping \
	-R2 TMP/dest.mapping > TMP/liftover.chain

test -s TMP/liftover.chain

wget -O - "${row.vcf}" |\
	java -jar \${JVARKIT_DIST}/vcfsetdict.jar -R TMP/ref.dict  --onNotFound SKIP |\
	bcftools annotate -O u --set-id +'%CHROM:%POS:%REF:%FIRST_ALT'  -x 'INFO,FILTER,QUAL' |\
	bcftools annotate -O u --set-id  '${row.name}:%ID' |\
	bcftools sort --max-mem "${task.memory.giga}G"  -T TMP -o TMP/${row.name}.bcf -O b



bcftools index TMP/${row.name}.bcf

mv TMP/${row.name}.bcf ./
mv TMP/${row.name}.bcf.csi ./
mv TMP/liftover.chain "${row.name}.chain"

bcftools view ${row.name}.bcf |\
	java -jar ${JVARKIT_DIST}/vcf2intervals.jar --bed --distance '50mb' --min-distance '250' |\
	awk -F '\t' -v P=\${PWD} '{printf("%s:%s-%s\t${row.name}\t%s/${row.name}.bcf\t%s/${row.name}.chain\t${row.reference}\\n",\$1,int(\$2)+1,\$3,P,P);}' > "${row.name}.bed"
"""
}

process APPLY_LIFTOVER {
tag "${name} ${interval}"
afterScript "rm -rf TMP"
memory "5g"
input:
	val(meta)
	tuple val(interval),val(name),val(vcf),val(chain),val(reference)
output:
	tuple val(name),path("lifted.bcf"),emit:output
script:
"""
hostname 1>&2
${moduleLoad("samtools jvarkit bcftools picard")}
mkdir -p TMP

bcftools view -o TMP/jeter1.vcf.gz -O z "${vcf}" "${interval}"
bcftools index TMP/jeter1.vcf.gz

java  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${PICARD_JAR} LiftoverVcf \
	--CHAIN "${chain}" \
	--INPUT TMP/jeter1.vcf.gz \
	--OUTPUT TMP/jeter2.vcf.gz \
	--REFERENCE_SEQUENCE "${reference}" \
	--REJECT TMP/reject.vcf.gz \
	--DISABLE_SORT true \
	--LOG_FAILED_INTERVALS false \
	--WRITE_ORIGINAL_ALLELES false \
	--WRITE_ORIGINAL_POSITION false 

bcftools annotate -x 'INFO,FILTER' -O u TMP/jeter2.vcf.gz |\
	bcftools sort --max-mem ${task.memory.giga}G -T TMP -O b -o "lifted.bcf"
bcftools index "lifted.bcf"
"""
}

process CONCATENATE {
input:
	val(meta)
	tuple val(name),val(vcfs)
output:
	tuple val(name),path("${name}.vcf.gz"),emit:output
script:
"""

"""
}
