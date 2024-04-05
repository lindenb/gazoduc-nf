include {moduleLoad;slurpJsonFile} from '../../modules/utils/functions.nf'


workflow {
	BUILD_DBSNP(file(params.json))
}

workflow BUILD_DBSNP {
	take:
		json
	main:
		json = slurpJsonFile(json)
		ch1_ch = Channel.of(json.sources).flatten()
		
		ch2_ch = DOWNLOAD_DICT(ch1_ch.
			filter{it.containsKey("fasta")}.
				map{[it.name,it.fasta]}
			)
		

		ch3_ch = DOWNLOAD_CHAIN(ch1_ch.
			filter{it.containsKey("chain")}.
				map{[it.name,it.chain]}
			)
		ch4_ch = DOWNLOAD_VCF(ch1_ch.flatMap{T->T.vcfs.collect{V->[T.name,V]}})	
		//VCF2INTERVALS(ch4_ch.output)
	}

process DOWNLOAD_DICT {
tag "${name} ${fasta}"
afterScript "rm -rf TMP"
input:
	tuple val(name),val(fasta)
output:
	tuple val(name),path("${name}.dict"),emit:output
script:
"""
hostname 1>&2
${moduleLoad("samtools")}
mkdir -p TMP
set -o pipefail

wget -q -O - "${fasta}" |\
	gunzip -c |\
	samtools dict -o TMP/ref.dict

mv  -v TMP/ref.dict "${name}.dict"
"""
}


process DOWNLOAD_CHAIN {
tag "${name} ${chain}"
afterScript "rm -rf TMP"
input:
	tuple val(name),val(chain)
output:
	tuple val(name),path("${name}.chain"),emit:output
script:
"""
hostname 1>&2
mkdir -p TMP
set -o pipefail

wget -q -O - "${chain}"  |\
	gunzip -c > TMP/liftover.chain

mv -v TMP/liftover.chain "${name}.chain"
"""
}

process DOWNLOAD_VCF {
tag "${name} ${vcf}"
afterScript "rm -rf TMP"
input:
	tuple val(name),val(vcf)
output:
	tuple val(name),path("${name}.bcf"),path("${name}.bcf.csi"),emit:output
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP
set -x
set -o pipefail

#
# we use vcf and NOT BCF because when there is not 'contig=' in the
# header, it cannot be converted to BCF
#
wget -q -O - "${vcf}"  |\\
	gunzip -c |\\
	cut -f1-8 |\\
	awk -F '\t' '/^##(INFO|FILTER|FORMAT)=/ {next;} /^#CHROM/ {printf("##source=${vcf}\\n");print;next;} /^#/ {print;next;} {ID=\$3; if(!(ID ~ /^[Rr][Ss][0-9]+\$/ )) {if(ID==".") {C=\$1;if(C ~ /^chr/) C=substr(C,4) ; ID=sprintf("%s_%s_%s",C,\$2,\$4);} ID=sprintf("${name}_%s",ID);} printf("%s\t%s\t%s\t%s\t%s\t.\t.\t.\\n",\$1,\$2,ID,\$4,\$5);}' |\
	bcftools view -O z -o TMP/jeter1.vcf.gz


#
# insert contigs if needed
#
if bcftools view --header-only TMP/jeter1.vcf.gz | grep '^##contig=' -m1
then
	bcftools view -O b -o TMP/jeter1.bcf TMP/jeter1.vcf.gz
else

	## try to build a dictionnary/fai , assume sorted on chromosome
	bcftools query -f '%CHROM\t%END\\n' TMP/jeter1.vcf.gz |\
		awk -F '\t' 'BEGIN{PREV_C="";PREV_LEN=1;} {if(PREV_C!=\$1 && PREV_C!="") {printf("%s\t%s\t%d\t1\t2\\n",PREV_C,PREV_LEN,NR);} PREV_C=\$1; PREV_LEN=\$2; } END {if(PREV_C!="") printf("%s\t%s\t%d\t1\t2\\n",PREV_C,PREV_LEN,NR);}' > TMP/jeter.fai


	bcftools reheader --fai TMP/jeter.fai -o TMP/jeter1.bcf --temp-prefix TMP/tmp TMP/jeter1.vcf.gz
fi


bcftools index -f TMP/jeter1.bcf

mv -v TMP/jeter1.bcf "${name}.bcf"
mv -v TMP/jeter1.bcf.csi "${name}.bcf.csi"
"""
}
/*
provess VCF2INTERVALS {
tag "${name} ${vcf.name}"
input:
	tuple val(name),path(vcf),path(csi)
output:
	tuple val(name),path("${name}.intervals"),path(vcf),path(csi),emit:output
script:
"""
hostname 1>&2
${moduleLoad("bcftools jvarkit")}
mkdir -p TMP
set -x

bcftools view "${vcf}" |\
        java -jar \${JVARKIT_DIST}/vcf2intervals.jar --bed --distance '50mb' --min-distance '250' |\
        awk -F '\t' '{printf("%s:%s-%s\\n",\$1,\int(\$2),\$3);}' > ${name}.intervals
"""
}*/

/*
process LIFTOVER_INTERVAL {
tag "${name} ${vcf.name} ${interval}"
memory '5g'
input:
	tuple val(name),val(interval),path(chain),path(vcf),path(csi)
output:
	path("lifted.bcf"),emit:output
script:
"""
hostname 1>&2
${moduleLoad("bcftools jvarkit")}
mkdir  -p TMP

bcftools view -O z -o TMP/jeter.vcf.gz --regions "${intervals}" "${bcf}" 
bcftools index -f -t TMP/jeter.vcf.gz

// build awk script to change chain
cat << EOF > TMP/jeter.awk
BEGIN	{
	OFS=" "
	IN=0;
EOF

bcftools query -l '${bcf}' |\\
	awk '{printf("dict1[\"%s\"]=\"%s\";\\n",\$1,\$1); if(\$1 ~ /^chr/) {C=substr(\$1,4);} else {C=sprintf("chr%s",\$1);} printf("if(!(\"%s\" in dict1)) {dict1[\"%s\"]=\"%s\";}\\n",C,\$1); }' >> TMP/jeter.awk
cut -f 1 "${reference}.fai" |\\
	awk '{printf("dict2[\"%s\"]=\"%s\";\\n",\$1,\$1); if(\$1 ~ /^chr/) {C=substr(\$1,4);} else {C=sprintf("chr%s",\$1);} printf("if(!(\"%s\" in dict2)) {dict2[\"%s\"]=\"%s\";}\\n",C,C,\$1); }' >> TMP/jeter.awk


cat << EOF >> TMP/jeter.awk
	}


/^chain/ {
	if((\$3 in dict1 && \$8 in dict2)) {
		IN=1;
		\$3 = dict1[\$3];
		\$8 = dict2[\$8];
		print;
		}
	else
		{
		IN=0;
		}
	next;
	}

	{
	if(IN==1) print;
	}
EOF

awk -f TMP/jeter.awk '${chain}' > TMP/jeter.chain

test -s TMP/jeter.chain

java  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${PICARD_JAR} LiftoverVcf \
	--CHAIN "TMP/jeter.chain" \
	--INPUT TMP/jeter.vcf.gz \
	--OUTPUT TMP/jeter2.vcf.gz \
	--REFERENCE_SEQUENCE "${reference}" \
	--REJECT TMP/reject.vcf.gz \
	--DISABLE_SORT true \
	--LOG_FAILED_INTERVALS false \
	--WRITE_ORIGINAL_ALLELES false \
	--WRITE_ORIGINAL_POSITION false 

bcftools annotate -x 'INFO,FILTER' -O u TMP/jeter2.vcf.gz |\
	bcftools sort --max-mem ${task.memory.giga}G -T TMP/x -O b -o "TMP/lifted.bcf"

bcftools index "TMP/lifted.bcf"

mv -v TMP/lifted.bcf ./
mv -v TMP/lifted.bcf.csi ./
"""
}
*/


/*
		// download fasta
		rsrc_ch.
		
			
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
*/

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
