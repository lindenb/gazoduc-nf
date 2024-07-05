workflow {

	beds_ch = Channel.fromPath(params.beds).
		splitText().
		map{file(it.trim())}

	bams_ch = Channel.fromPath(params.bams).
		splitText().
		map{it.trim()}.
		flatMap{[it, it.endsWith(".cram")?it+".crai":it+".bai"]}.
		map{file(it,checkIfExists:true)}.
		toSortedList()

	genome_ch = [file(params.fasta) , file(params.fasta+".fai"), file(""+file(params.fasta).getParent()+"/"+file(params.fasta).getBaseName()+".dict" ) ]	
	if(params.dbsnp.equals("NO_FILE"))
		{
		dbsnp_ch = [file("NO_DBSNP"), file("NO_DBSNP_TBI")]
		}
	else
		{
		dbsnp_ch = [file(params.dbsnp), file(params.dbsnp+".tbi")]
		}
	hc_ch = HC_BY_BED( genome_ch, dbsnp_ch, bams_ch, beds_ch)
	VCF_CONCAT(hc_ch.flatten().collect())
	}

process HC_BY_BED {
tag "${bed.name}"
label "process_low"
memory "5g"
afterScript "rm -rf TMP"
cpus 4
input:
	tuple path(fasta),path(fai),path(dict)
	tuple path(dbsnp),path(dbsnp_tbi)
	path("BAMS/*")
	path(bed)
output:
	tuple path("${bed.getBaseName()}.vcf.gz"),path("${bed.getBaseName()}.vcf.gz.tbi"),emit:output
script:
	def min_file_split=25
"""
hostname 1>&2
module load gatk/0.0.0 bcftools
mkdir -p TMP
set -x

find ./BAMS -name "*.bam" -o -name "*.cram" | sort > TMP/bams.list

i=1
cat TMP/bams.list | while read B
do

   gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" HaplotypeCaller \\
     -L "${bed}" \\
     -R "${fasta}" \\
     -I "\${B}" \\
     -ERC GVCF \\
     ${dbsnp.name.equals("NO_DBSNP")?"":"--dbsnp ${dbsnp}"} \\
     ${(params.mapq as Integer)<1?"":" --minimum-mapping-quality "+params.mapq} \\
     -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \\
     -O "TMP/jeter.\${i}.g.vcf.gz"

     echo "TMP/jeter.\${i}.g.vcf.gz" >> TMP/gvcfs.list
     i=\$((i+1))		

done

SQRT=`awk 'END{X=NR;if(X <= ${min_file_split}){print(X);} else {z=sqrt(X); print (z==int(z)?z:int(z)+1);}}' TMP/gvcfs.list`

split -a 9 --additional-suffix=.list --lines=\${SQRT} TMP/gvcfs.list  TMP/split.

i=1
find ./TMP -type f -name "split*.list" | while read F
do
	gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" \\
		CombineGVCFs \\
		-R "${fasta}" \\
		-L "${bed}" \\
		-V "\${F}" \\
		-O "TMP/combined.\${i}.g.vcf.gz" \\
		-G StandardAnnotation \\
		-G AS_StandardAnnotation

	 echo "TMP/combined.\${i}.g.vcf.gz" >> TMP/gvcfs.list
	 i=\$((i+1))
done


	gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" \\
		CombineGVCFs \\
		-R "${fasta}" \\
		-L "${bed}" \\
		-V "TMP/gvcfs.list" \\
		-O TMP/combined.all.g.vcf.gz \\
		-G StandardAnnotation \\
		-G AS_StandardAnnotation

	gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" \\
		GenotypeGVCFs \\
		-R "${fasta}" \\
		-L "${bed}" \\
		-O TMP/genotyped.vcf.gz \\
		-V TMP/combined.all.g.vcf.gz \\
		-G StandardAnnotation \\
		-G AS_StandardAnnotation



mv -v TMP/genotyped.vcf.gz "${bed.getBaseName()}.vcf.gz"
mv -v TMP/genotyped.vcf.gz.tbi "${bed.getBaseName()}.vcf.gz.tbi"
"""
}

process VCF_CONCAT {
label "process_low"
memory "5g"
input:
	path("VCFS/*")
output:
	tuple path("output.bcf"),path("output.bcf.csi"),emit:output
script:
	
"""
hostname 1>&2
module load gatk/0.0.0 bcftools
mkdir -p TMP
	
find ./VCFS -name "*.vcf.gz" > TMP/jeter.list

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" \\
		MergeVcfs \\
		I=TMP/jeter.list \\
		O=TMP/output.vcf.gz

bcftools view -O b9 -o TMP/output.bcf TMP/output.vcf.gz
bcftools index -f TMP/output.bcf
mv TMP/output.bcf ./
mv TMP/output.bcf.csi ./
"""
}
