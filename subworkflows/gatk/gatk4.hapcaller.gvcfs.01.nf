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


include {parseBoolean;isBlank;moduleLoad;getVersionCmd} from './../../modules/utils/functions.nf'
include {gatkGetArgumentsForCombineGVCFs;gatkGetArgumentsForGenotypeGVCF} from './gatk.hc.utils.nf'
include {SAMTOOLS_SAMPLES} from '../samtools/samtools.samples.03.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {GATK_PARALLEL_COMBINE_GVCFS} from './gatk.parallel.combine.gvcfs.01.nf'

workflow GATK4_HAPCALLER_GVCFS_01 {
        take:
		meta
		genomeId
                bams
                beds
		pedigree
        main:
                version_ch = Channel.empty()


		samples_bams_ch = SAMTOOLS_SAMPLES(
			[allow_multiple_references:true,allow_duplicate_samples:true],
			bams
			)
		version_ch = version_ch.mix(samples_bams_ch.version)


		each_bed = beds.splitText().map{it.trim()}

		//each_bam.dump()		

		hc_input_ch = each_bed.combine(samples_bams_ch.rows).map{T->[
			"bed":T[0],
			"sample":T[1].sample,
			"bam":T[1].bam,
			"new_sample":T[1].new_sample,
			"genomeId": T[1].genomeId,
			"pedigree": (pedigree.name.equals("NO_FILE")?"":pedigree.toRealPath()),
			]}
		hapcaller_ch = HAPCALLER_GVCF([:],genomeId,hc_input_ch)
		version_ch = version_ch.mix(hapcaller_ch.version.first())

		by_bed_ch = hapcaller_ch.bedvcf.map{T->[T[0].bed,T[1]]}.groupTuple()

		split_gvcfs_ch = SPLIT_GVCF_BLOCKS_PER_CONTIG(by_bed_ch)
		version_ch = version_ch.mix(split_gvcfs_ch.version.first())


		find_gvcfs_ch = FIND_GVCF_BLOCKS([:], genomeId, split_gvcfs_ch.output.splitCsv(header:true,sep:'\t',strip:true))
		version_ch = version_ch.mix(find_gvcfs_ch.version.first())

		each_chunk = find_gvcfs_ch.output.splitCsv(header:true,sep:'\t',strip:true)

		if( parseBoolean(params.parallel_combine_gvcf) )  {
			genotyped_ch = GATK_PARALLEL_COMBINE_GVCFS([:], genomeId, pedigree, each_chunk)
			version_ch = version_ch.mix(genotyped_ch.version.first())
			}
		else {
			genotyped_ch = GENOTYPE_GVCFS_02([:], genomeId, pedigree, each_chunk)
			version_ch = version_ch.mix(genotyped_ch.version.first())
			}

		version_ch = MERGE_VERSION("gatk4",version_ch.collect())
	emit:
		version = version_ch
		region_vcf = genotyped_ch.region_vcf
}



process HAPCALLER_GVCF {
tag "bed:${file(row.bed).name} genomeId:${row.genomeId} bam:${file(row.bam).name}"
cache 'lenient'
memory {task.attempt==1?'10G':'20G'}
errorStrategy 'retry'
maxRetries 5
cpus 1
afterScript 'rm -rf TMP'
input:
	val(meta)
	val(genomeId)
	val(row)
output:
        tuple val(row),path("hapcaller.g.vcf.gz"),emit:bedvcf
        path("hapcaller.g.vcf.gz.tbi"),emit:index
        path("version.xml"),emit:version
script:
	def bam = row.bam?:""
	def sample = row.sample?:""
	def new_sample = row.new_sample?:""
	def extraHC = params.gatk.haplotypecaller.args?:""
	def pedigree = row.pedigree?:""
	def mapq = params.gatk.haplotypecaller.mapq?:-1
	def bams = row.bams?:""
	def reference = params.genomes[genomeId].fasta
	def bam_reference = params.genomes[row.genomeId].fasta
	def bed = row.bed?:""
	def dbsnp =  params.genomes[row.genomeId].dbsnp?:""
	def fix_tp_ap =  params.gatk.haplotypecaller.fix_tp_ap?:false
"""
hostname 1>&2
${moduleLoad("gatk4 bcftools jvarkit samtools")}

   set -o pipefail
   set -x
   mkdir -p TMP

   if [ "${genomeId}" != "${row.genomeId}" ] ; then
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/bedrenamechr.jar -R "${bam_reference}" --column 1 --convert SKIP "${bed}" > TMP/fixed.bed
	if [ ! -s TMP/fixed.bed ] ; then
		tail -n 1 "\${REF}.fai" | awk -F '\t' '{printf("%s\t0\t1\\n",\$1);}' > TMP/fixed.bed
   	fi 
   else
	ln -s "${bed}" TMP/fixed.bed
   fi

   # 20231116 https://gatk.broadinstitute.org/hc/en-us/community/posts/11440622639387-Unable-to-trim-uncertain-bases-without-flow-order-information
   #
   if ${fix_tp_ap} ; then
	samtools view -M -L TMP/fixed.bed -T '${bam_reference}' -h '${bam}' |\
		sed 's/\ttp\\:A\\:P\t/\t/' |\
		samtools view -O BAM -o TMP/jeter.bam
	samtools index TMP/jeter.bam
		
   fi 

   gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" HaplotypeCaller \
     -I "${fix_tp_ap?"TMP/jeter.bam":bam}" \
     -ERC GVCF \
     --seconds-between-progress-updates 600 \
     ${!isBlank(dbsnp) && genomeId.equals(row.genomeId) ?"--dbsnp ${dbsnp}":""}   \
     -L TMP/fixed.bed \
     -R "${bam_reference}" \
     ${isBlank(pedigree)?"":" --pedigree '${pedigree}'"} \
     ${(mapq as Integer)<1?"":" --minimum-mapping-quality '${mapq}'"} \
     -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
     ${extraHC} \
     -O "TMP/jeter.g.vcf.gz"

   # reference is not the main reference
   if [ "${reference}" != "${bam_reference}" ] ; then

	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/vcfsetdict.jar -n SKIP -R "${reference}" TMP/jeter.g.vcf.gz > TMP/jeter2.vcf
	
	bcftools sort -T TMP  --max-mem "${task.memory.giga}G" -O z -o TMP/jeter.g.vcf.gz TMP/jeter2.vcf
	bcftools index --force --tbi TMP/jeter.g.vcf.gz
	rm TMP/jeter2.vcf

	# annotate dbsnp
	if ${!isBlank(dbsnp)} ] ; then
		# bcftools annotate is just too slow

		gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"  VariantAnnotator \
			-R "${reference}" \
			--dbsnp "${dbsnp}" \
			-V TMP/jeter.g.vcf.gz \
			-O TMP/jeter2.g.vcf.gz
		mv TMP/jeter2.g.vcf.gz TMP/jeter.g.vcf.gz
		mv TMP/jeter2.g.vcf.gz.tbi TMP/jeter.g.vcf.gz.tbi
	fi

   fi
 
   # must rename sample
   if [ "${sample}" != "${new_sample}" ] ; then

	gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"  RenameSampleInVcf \
		-INPUT TMP/jeter.g.vcf.gz \
		-OUTPUT TMP/jeter2.g.vcf.gz \
		-NEW_SAMPLE_NAME "${new_sample}"

	bcftools index --force --tbi TMP/jeter2.g.vcf.gz
	mv TMP/jeter2.g.vcf.gz TMP/jeter.g.vcf.gz
	mv TMP/jeter2.g.vcf.gz.tbi TMP/jeter.g.vcf.gz.tbi
   fi



   mv TMP/jeter.g.vcf.gz  hapcaller.g.vcf.gz
   mv TMP/jeter.g.vcf.gz.tbi  hapcaller.g.vcf.gz.tbi


##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">haplotype caller in GVCF mode</entry>
        <entry key="bam">${bam}</entry>
        <entry key="sample">${sample}</entry>
        <entry key="new.sample">${new_sample}</entry>
        <entry key="reference">${reference}</entry>
        <entry key="bam_reference">${bam_reference}</entry>
        <entry key="dbsnp">${dbsnp}</entry>
        <entry key="pedigree">${pedigree}</entry>
	<entry key="gatk.version">\$( gatk HaplotypCaller --version 2> /dev/null  | paste -s -d ' ')</entry>
	<entry key="gatk.cmd"><code>\$(bcftools view --header-only hapcaller.g.vcf.gz | grep "^##GATKCommandLine=<ID=HaplotypeCaller" -m1 | tr "<>&" "_" | cat)</code></entry>
</properties>
EOF
"""
}


process SPLIT_GVCF_BLOCKS_PER_CONTIG {
executor "local"
tag "${file(bed).name} N=${L.size()}"
memory {task.attempt==1?"5g":(task.attempt==2?"15g":"30g")}
errorStrategy "retry"
maxRetries 3
afterScript "rm -r jeter.list"
input:
	tuple val(bed),val(L)
output:
	path("contig.gvcfs.bed.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def sqrt = (L.size() < 100 ? L.size() : Math.max(1,(int)Math.sqrt(L.size())))
"""
hostname 1>&2
set -o pipefail

mkdir -p GVCFS

cat << EOF > gvcfs.list
${L.join("\n")}
EOF

split -a 9 --additional-suffix=.list --lines=${sqrt} gvcfs.list  GVCFS/cluster0.


find \${PWD}/GVCFS -type f -name "cluster0.*.list" > split.gvcfs.txt
test -s split.gvcfs.txt

cut -f 1 "${bed}" | sort -T . | uniq |\
	awk -v P=\${PWD} 'BEGIN{printf("contig\tgvcf_list\tgvcf_split\tbed\\n");} {printf("%s\t%s/gvcfs.list\t%s/split.gvcfs.txt\t${bed}\\n",\$1,P,P);}' > contig.gvcfs.bed.tsv


# prevent timestamp problem
sleep 10

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">make a list of all GVCFS for a bed</entry>
        <entry key="bed">${bed}</entry>
        <entry key="sqrt">${sqrt}</entry>
        <entry key="N(gvcfs)">${L.size()}</entry>
</properties>
EOF
"""

}



process FIND_GVCF_BLOCKS {
tag "${row.contig}"
memory "10g"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(genomeId)
	val(row)
output:
	path("bed_samplemap.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def contig = row.contig
	def blocksize= params.blocksize?:"1mb"
	def mergesize = params.mergesize?:"100"
	def reference = params.genomes[genomeId].fasta

if(parseBoolean(params.use_whole_block))
"""
hostname 1>&2
${moduleLoad("jvarkit")}

mkdir -p TMP

awk -F '\t' 'BEGIN{printf("interval\tgvcf_split\tbed\\n");X=-1;Y=-1;} (\$1=="${contig}") {B=int(\$2);E=int(\$3);if(X==-1) {X=B;Y=E;} if(B<X) X=B; if(E>Y) Y=E;} END {printf("${contig}:%d-%d\t${row.gvcf_split}\t${row.bed}\\n",X+1,Y);}' '${row.bed}' > bed_samplemap.tsv

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">use whome blocks</entry>
        <entry key="contig">${contig}</entry>
</properties>
EOF
"""
else
"""
hostname 1>&2
${moduleLoad("jvarkit")}

mkdir -p TMP

java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/findgvcfsblocks.jar \
	-T TMP \
	--bed "${row.bed}" \
	--contig "${contig}" \
	--block-size "${blocksize}" \
	--merge-size "${mergesize}" \
	-o "TMP/jeter.interval_list" \
	"${row.gvcf_list}"

awk -F '\t' 'BEGIN{printf("interval\tgvcf_split\tbed\\n");} /^@/ {next;} {printf("%s:%d-%s\t${row.gvcf_split}\t${row.bed}\\n",\$1,\$2,\$3);}' TMP/jeter.interval_list > bed_samplemap.tsv

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">find gvcfs blocks using findgvcfsblocks</entry>
        <entry key="contig">${contig}</entry>
        <entry key="block.size">${blocksize}</entry>
        <entry key="N(output)">\$(wc -l < bed_samplemap.tsv)</entry>
</properties>
EOF
"""
}



process GENOTYPE_GVCFS_02 {
tag "${row.interval}"
memory {task.attempt <2 ? "15g":"60g"}
memory "10g"
cpus "3"
errorStrategy  'retry'
maxRetries 3
afterScript 'rm -rf  TMP BEDS jeter* database tmp_read_resource_*.config'
input:
	val(meta)
	val(genomeId)
	path(pedigree)
	val(row)
output:
        tuple val("${row.interval}"),path("genotyped.bcf"),emit:region_vcf
        path("genotyped.bcf.csi"),emit:index
	path("version.xml"),emit:version
script:
	 def reference = params.genomes[genomeId].fasta
	 def otherOpts1 =  gatkGetArgumentsForCombineGVCFs(params.plus("genomeId":genomeId))
	 def otherOpts2 =  gatkGetArgumentsForGenotypeGVCF(params.plus("genomeId":genomeId,"pedigree":pedigree))
         def region = row.interval
"""
hostname 1>&2
${moduleLoad("gatk4 bcftools")}

touch vcfs.list
mkdir BEDS TMP

j=1
cat "${row.gvcf_split}" | while read V
do

	gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" CombineGVCFs \
		${otherOpts1} \
		-L "${region}" \
		-V "\${V}"  \
		-O "TMP/combine0.\${j}.g.vcf.gz"

	echo "TMP/combine0.\${j}.g.vcf.gz" >> TMP/combine0.list
	j=\$((j+1))
	
done

test -s TMP/combine0.list

if [ "\${j}" == "2" ] ; then

	mv TMP/combine0.1.g.vcf.gz TMP/combine1.g.vcf.gz
	mv TMP/combine0.1.g.vcf.gz.tbi TMP/combine1.g.vcf.gz.tbi

else

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" CombineGVCFs \
        ${otherOpts1} \
        -L "${region}" \
	-V "TMP/combine0.list"  \
	-O "TMP/combine1.g.vcf.gz"
fi

rm -f TMP/combine0.list


gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" GenotypeGVCFs  \
      ${otherOpts2} \
      -L "${region}" \
      -V TMP/combine1.g.vcf.gz \
      --seconds-between-progress-updates 60 \
      -O "TMP/jeter.vcf.gz"


bcftools view  --threads ${task.cpus}  -O b -o genotyped.bcf TMP/jeter.vcf.gz
bcftools index  --threads ${task.cpus}  genotyped.bcf

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">combine and call gvcfs</entry>
        <entry key="region">${region}</entry>
	<entry key="gatk.version">\$( gatk HaplotypCaller --version 2> /dev/null  | paste -s -d ' ')</entry>
</properties>
EOF
"""
}



