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
		

		ch3_ch = DOWNLOAD_CHAIN(
			json.name,
			ch1_ch.map{it.name}
			)
		ch4_ch = DOWNLOAD_VCF(ch1_ch.map{T->{
			if(T.containsKey("vcf") && !T.containsKey("vcfs")) {
				return T.plus("vcfs":[T.vcf]);
				}
			else
				{
				return T;
				}
			}}.flatMap{T->T.vcfs.collect{V->[T.name,V]}})	

		ch4_ch.output.view()
		ch5_ch = VCF2INTERVALS(ch4_ch.output)
		
		ch6_ch = ch5_ch.output.splitText().
			map{T->{T[1]=T[1].trim(); return T;}}.
			join(ch3_ch.output)
		
		ch7_ch = LIFTOVER_INTERVAL(json.name,json.fasta,ch6_ch)
		CONCATENATE(json.name, ch7_ch.output.flatten().collect())
			
	}


process DOWNLOAD_CHAIN {
tag "${name}->${hgTo}"
afterScript "rm -rf TMP"
input:
	val(hgTo)
	val(name)
output:
	tuple val(name),path("${name}.chain"),optional:true,emit:output
script:        
	def chainName = "${name}To${hgTo.substring(0,1).toUpperCase()}${hgTo.substring(1)}.over.chain"
	def url = "http://hgdownload.cse.ucsc.edu/goldenpath/${name}/liftOver/${chainName}.gz"
"""
hostname 1>&2
mkdir -p TMP
set -o pipefail

if wget --spider "${url}"
then

	wget -q -O - "${url}"  |\
		gunzip -c > TMP/liftover.chain

	test -s TMP/liftover.chain
	mv -v TMP/liftover.chain "${name}.chain"
fi
"""
}

process DOWNLOAD_VCF {
tag "${name} ${vcf}"
afterScript "rm -rf TMP"
input:
	tuple val(name),val(vcf)
output:
	tuple val(name),path("${name}.vcf.gz"),path("${name}.vcf.gz.tbi"),emit:output
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP
set -x

#
# we use vcf and NOT BCF because when there is not 'contig=' in the
# header, it cannot be converted to BCF
#
wget -q -O - "${vcf}"  |\\
	bcftools view |\\
	awk -F '\t' '/^##(INFO|FILTER|FORMAT)=/ {next;} /^#CHROM/ {printf("##source.${name}=${vcf}\\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\\n");N=0;next;} /^#/ {print;next;} {${(params.limit as int)> 0 ?"if(N> ${params.limit}) {exit 0;};N++;":""}printf("%s\t%s\t.\t%s\t%s\t.\t.\t.\\n",\$1,\$2,\$4,\$5);}' |\\
	bcftools view -O z -o TMP/jeter1.vcf.gz


#
# insert contigs if needed
#
if bcftools view --header-only TMP/jeter1.vcf.gz | grep '^##contig=' -m1
then
	mv TMP/jeter1.vcf.gz TMP/jeter2.vcf.gz
else

	## try to build a dictionnary/fai , assume sorted on chromosome
	bcftools query -f '%CHROM\\n' TMP/jeter1.vcf.gz |\
		uniq |\\
		awk -F '\t' '{printf("%s\t1000000000\t%d\t1\t2\\n",\$1,NR);}' > TMP/jeter.fai


	bcftools reheader --fai TMP/jeter.fai -o TMP/jeter2.vcf.gz --temp-prefix TMP/tmp TMP/jeter1.vcf.gz
fi

mv TMP/jeter2.vcf.gz TMP/jeter1.vcf.gz
bcftools index -t -f TMP/jeter1.vcf.gz

mv -v TMP/jeter1.vcf.gz "${name}.vcf.gz"
mv -v TMP/jeter1.vcf.gz.tbi "${name}.vcf.gz.tbi"
"""
}

process VCF2INTERVALS {
tag "${name} ${vcf.name}"
afterScript "rm -rf TMP"
memory '1G'
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

bcftools view "${vcf}" |\\
        java -jar -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP \${JVARKIT_DIST}/vcf2intervals.jar --bed --distance '50mb' --min-distance '250' |\\
        awk -F '\t' '{printf("%s:%s-%s\\n",\$1,int(\$2)+1,\$3);}' > TMP/jeter.intervals

mv -v TMP/jeter.intervals ${name}.intervals
"""
}


process LIFTOVER_INTERVAL {
tag "${name}  ${interval}"
afterScript "rm -rf TMP"
memory '5g'
input:
	val(toName)
	val(fasta)
	tuple val(name),val(interval),path(vcf),path(csi),path(chain)
output:
	path("${name}.${interval.replaceAll("[^A-Za-z0-9]+","_")}.lifted.*"),emit:output
script:
"""
hostname 1>&2
${moduleLoad("bcftools jvarkit picard")}
mkdir  -p TMP
set -x
test -s "${fasta}.fai"

bcftools view -O z -o TMP/jeter.vcf.gz --regions "${interval}" "${vcf}" 
bcftools index -f -t TMP/jeter.vcf.gz


#
# build awk script to change chain
#
cat << EOF > TMP/jeter.awk
BEGIN	{
	OFS=" "
	IN=0;
EOF


bcftools index -s '${vcf}' |\\
	awk -F '\t' '{printf("\tdict1[\\\"%s\\\"]=\\\"%s\\\";\\n",\$1,\$1); if(\$1 ~ /^chr/) {C=substr(\$1,4);} else {C=sprintf("chr%s",\$1);} printf("\tif(!(\\\"%s\\\" in dict1)) {dict1[\\\"%s\\\"]=\\\"%s\\\";}\\n",C,C,\$1); }' >> TMP/jeter.awk

cut -f 1 "${fasta}.fai" |\\
	awk '{printf("\tdict2[\\\"%s\\\"]=\\\"%s\\\";\\n",\$1,\$1); if(\$1 ~ /^chr/) {C=substr(\$1,4);} else {C=sprintf("chr%s",\$1);} printf("\tif(!(\\\"%s\\\" in dict2)) {dict2[\\\"%s\\\"]=\\\"%s\\\";}\\n",C,C,\$1); }' >> TMP/jeter.awk


cat << 'EOF' >> TMP/jeter.awk
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


java  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${PICARD_JAR} LiftoverVcf \\
	--CHAIN "TMP/jeter.chain" \\
	--INPUT TMP/jeter.vcf.gz \\
	--OUTPUT TMP/jeter2.vcf.gz \\
	--REFERENCE_SEQUENCE "${fasta}" \\
	--REJECT TMP/reject.vcf.gz \\
	--DISABLE_SORT true \\
	--LOG_FAILED_INTERVALS false \\
	--WRITE_ORIGINAL_ALLELES false \\
	--WRITE_ORIGINAL_POSITION true 

bcftools annotate --set-id '${name}_%INFO/OriginalContig\\_%INFO/OriginalStart'  -O u TMP/jeter2.vcf.gz |\\
	bcftools annotate -x 'INFO,FILTER' -O u |\
	bcftools sort --max-mem ${task.memory.giga}G -T TMP/x -O b -o "TMP/lifted.bcf"

bcftools index "TMP/lifted.bcf"

mv -v TMP/lifted.bcf "${name}.${interval.replaceAll("[^A-Za-z0-9]+","_")}.lifted.bcf"
mv -v TMP/lifted.bcf.csi "${name}.${interval.replaceAll("[^A-Za-z0-9]+","_")}.lifted.bcf.csi"
"""
}


process CONCATENATE {
afterScript "rm -rf TMP"
input:
	val(hgTo)
	path(vcfs)
output:
	path("${params.prefix}${hgTo}.vcf.gz"),emit:output
script:
"""

hostname 1>&2
${moduleLoad("bcftools jvarkit picard")}
mkdir  -p TMP

bcftools  concat --remove-duplicates -a -O b -o TMP/concat.vcf.gz *.bcf
bcftools index -f -t TMP/concat.vcf.gz

mv TMP/concat.vcf.gz "${params.prefix}${hgTo}.vcf.gz"
"""
}
