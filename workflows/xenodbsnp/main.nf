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
			ch1_ch.filter{it.liftover==true}.map{it.name}
			)
		ch4_ch = DOWNLOAD_VCF(ch1_ch.
			filter{it.liftover==true}.
			map{T->{
			if(T.containsKey("vcf") && !T.containsKey("vcfs")) {
				return T.plus("vcfs":[T.vcf]);
				}
			else
				{
				return T;
				}
			}}.flatMap{T->T.vcfs.collect{V->[T.name,T.priority,V]}})	

		ch5_ch = VCF2INTERVALS(ch4_ch.output)
		
		ch6_ch = ch5_ch.output.splitText().
			map{T->{T[2]=T[2].trim(); return T;}}.
			join(ch3_ch.output)
		
		ch7_ch = LIFTOVER_INTERVAL(json.name,json.fasta,ch6_ch)


		ch8_ch = DOWNLOAD_VCF_NO_LIFT(
			json.fasta,
			ch1_ch.
                        filter{it.liftover==false}.
			map{T->{
			if(T.containsKey("vcf") && !T.containsKey("vcfs")) {
				return T.plus("vcfs":[T.vcf]);
				}
			else
				{
				return T;
				}
			}}.flatMap{T->T.vcfs.collect{V->[T.name,T.priority,V]}}
			)
	

		CONCATENATE(json.name, json.fasta, ch7_ch.output.mix(ch8_ch.output).
			map{""+it[0]+"\t"+it[1]}.
			collect())
			
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
	tuple val(name),val(priority),val(vcf)
output:
	tuple val(name),val(priority),path("${name}.vcf.gz"),path("${name}.vcf.gz.tbi"),emit:output
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
	tuple val(name),val(priority),path(vcf),path(csi)
output:
	tuple val(name),val(priority),path("${name}.intervals"),path(vcf),path(csi),emit:output
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
	tuple val(name),val(priority),val(interval),path(vcf),path(csi),path(chain)
output:
	tuple val(priority),path("${name}.${interval.replaceAll("[^A-Za-z0-9]+","_")}.lifted.vcf.gz"),emit:output
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
	bcftools annotate -x 'INFO,FILTER' -O v |\
	awk -F '\t' '/^#/ {print;next} {OFS="\t";ID=\$3;N=split(ID,a,/_/);if(N==3 && a[2] ~ /^chr/) {ID=sprintf("%s_%s_%s",a[1],substr(a[2],4),a[3]);\$3=ID;}print;}' |\
	bcftools sort --max-mem ${task.memory.giga}G -T TMP/x -O z -o "TMP/lifted.vcf.gz"

bcftools index -f -t "TMP/lifted.vcf.gz"

mv -v TMP/lifted.vcf.gz "${name}.${interval.replaceAll("[^A-Za-z0-9]+","_")}.lifted.vcf.gz"
mv -v TMP/lifted.vcf.gz.tbi "${name}.${interval.replaceAll("[^A-Za-z0-9]+","_")}.lifted.vcf.gz.tbi"
"""
}


process DOWNLOAD_VCF_NO_LIFT {
tag "${name} ${vcf}"
afterScript "rm -rf TMP"
input:
	val(fasta)
	tuple val(name),val(priority),val(vcf)
output:
	tuple val(priority),path("${name}.vcf.gz"),emit:output
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP
set -x

cat << 'EOF' > TMP/jeter.awk
BEGIN	{
	FS="\t";
	OFS="\t";
EOF

test -f '${fasta}.fai'

awk -F '\t' '{printf("\tdict1[\\\"%s\\\"]=\\\"%s\\\";\\n",\$1,\$1); if(\$1 ~ /^chr/) {C=substr(\$1,4);} else {C=sprintf("chr%s",\$1);} printf("\tif(!(\\\"%s\\\" in dict1)) {dict1[\\\"%s\\\"]=\\\"%s\\\";}\\n",C,C,\$1); }' '${fasta}.fai'  >> TMP/jeter.awk

cat << 'EOF' >> TMP/jeter.awk
	}

/^##(INFO|FILTER|FORMAT)=/ {
	next;
	}

/^#CHROM/ {
	printf("##source.${name}=${vcf}\\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\\n");
	N=0;
	next;
	}

/^#/	{
	print;
	next;
	}

	{
	if(!(\$1 in dict1)) next;
	${(params.limit as int)> 0 ?"if(N> ${params.limit}) {exit 0;};N++;":""}
	C = dict1[\$1];
	printf("%s\t%s\t",C,\$2);
	ID= \$3;
	if(!(ID ~ /^[rR][sS][0-9]+\$/ )) {
		if(C ~ /^chr/) {
			ID= sprintf("${name}_%s_%s",substr(C,4),\$2);
			}
		else
			{
			ID= sprintf("${name}_%s_%s",C,\$2);
			}
		}
	printf("%s\t%s\t%s\t.\t.\t.\\n",ID,\$4,\$5);
	}
EOF


wget -O - "${vcf}" |\\
	gunzip -c |\\
	awk -f TMP/jeter.awk |\\
	bcftools view -O z -o TMP/jeter1.vcf.gz

bcftools index -t -f TMP/jeter1.vcf.gz

mv -v TMP/jeter1.vcf.gz "${name}.vcf.gz"
mv -v TMP/jeter1.vcf.gz.tbi "${name}.vcf.gz.tbi"
"""
}


process CONCATENATE {
tag "N=${vcfs.size()}"
afterScript "rm -rf TMP"
input:
	val(hgTo)
	val(fasta)
	val(vcfs)
output:
	path("${params.prefix}${hgTo}.vcf.gz"),emit:output
script:
"""

hostname 1>&2
${moduleLoad("bcftools jvarkit")}
mkdir  -p TMP/TMP

cat << EOF > TMP/jeter.list
${vcfs.join("\n")}
EOF

test -s TMP/jeter.list

java -version 1>&2                                  
                                                                                                               
cp -v "${moduleDir}/Minikit.java" TMP/TMP/Minikit.java                                                        
javac -d TMP/TMP -cp \${JVARKIT_DIST}/jvarkit.jar -sourcepath 'TMP/TMP' 'TMP/TMP/Minikit.java'                            
jar cvf TMP/minikit.jar -C TMP/TMP .                                                                           

java -cp \${JVARKIT_DIST}/jvarkit.jar:TMP/minikit.jar Minikit '${fasta}' TMP/jeter.list |\
	bcftools view -O z -o TMP/concat.vcf.gz

mv -v TMP/concat.vcf.gz "${params.prefix}${hgTo}.vcf.gz"
"""
}
