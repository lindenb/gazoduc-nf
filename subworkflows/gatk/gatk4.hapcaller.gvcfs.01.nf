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

def gazoduc = gazoduc.Gazoduc.getInstance()

gazoduc.make("parallel_combine_gvcf",false).
        description("Run a different process for each step of combineGVCFs + GenotypeGVCFs. Default is to CombineGVCF+GenotypeGVCS in the same process").
	setBoolean().
        put()


gazoduc.make("use_whole_block",false).
        description("Do NOT use jvarkit/findgvcfsblocks but use the whole interval. Useful when the BED is a set of small exons (WES).").
	setBoolean().
        put()


gazoduc.make("blocksize","1mb").
        description("argument for jvarkit/findgvcfsblocks block-size").
        put()

gazoduc.make("mergesize","100").
        description("argument for jvarkit/findgvcfsblocks merge-size").
        put()

gazoduc.make("extraHC","").
        description("extra arguments for haplotype caller").
        put()


include {parseBoolean;isBlank;moduleLoad;getVersionCmd} from './../../modules/utils/functions.nf'
include {gatkGetArgumentsForCombineGVCFs;gatkGetArgumentsForGenotypeGVCF} from './gatk.hc.utils.nf'
include {SAMTOOLS_SAMPLES_01} from '../samtools/samtools.samples.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {GATK_PARALLEL_COMBINE_GVCFS} from './gatk.parallel.combine.gvcfs.01.nf'

workflow GATK4_HAPCALLER_GVCFS_01 {
        take:
                meta
                reference
		references
                bams
                beds
		pedigree
        main:
                version_ch = Channel.empty()


		samples_bams_ch = SAMTOOLS_SAMPLES_01(meta.plus(["with_header":true,"allow_duplicate_samples":true,"allow_multiple_references":true]), reference, references, bams)
		version_ch = version_ch.mix(samples_bams_ch.version)


		each_bed = beds.splitText().map{it.trim()}
		each_bam = samples_bams_ch.output.splitCsv(header: true, sep : '\t')

		//each_bam.dump()		

		hc_input_ch = each_bed.combine(each_bam).map{T->[
			"bed":T[0],
			"sample":T[1].sample,
			"bam":T[1].bam,
			"new_sample":T[1].new_sample,
			"reference": reference,
			"dbsnp": (meta.dbsnp?:""),
			"extraHC": (meta.extraHC?:""),
			"pedigree": (pedigree.name.equals("NO_FILE")?"":pedigree.toRealPath()),
			"mapq": (meta.mapq?:"-1"),
			"bam_reference": T[1].reference
			]}
		hapcaller_ch = HAPCALLER_GVCF(meta, hc_input_ch)
		version_ch = version_ch.mix(hapcaller_ch.version.first())

		by_bed_ch = hapcaller_ch.bedvcf.map{T->[T[0].bed,T[1]]}.groupTuple()

		split_gvcfs_ch = SPLIT_GVCF_BLOCKS_PER_CONTIG(meta,by_bed_ch)
		version_ch = version_ch.mix(split_gvcfs_ch.version.first())


		find_gvcfs_ch = FIND_GVCF_BLOCKS(meta,split_gvcfs_ch.output.splitCsv(header:true,sep:'\t',strip:true))
		version_ch = version_ch.mix(find_gvcfs_ch.version.first())

		each_chunk = find_gvcfs_ch.output.splitCsv(header:true,sep:'\t',strip:true)

		if( parseBoolean(meta.parallel_combine_gvcf) )  {
			genotyped_ch = GATK_PARALLEL_COMBINE_GVCFS(meta,reference, pedigree, each_chunk)
			version_ch = version_ch.mix(genotyped_ch.version.first())
			}
		else {
			genotyped_ch = GENOTYPE_GVCFS_02(meta,reference, pedigree, each_chunk)
			version_ch = version_ch.mix(genotyped_ch.version.first())
			}

		version_ch = MERGE_VERSION(meta, "gatk4", "call bams using gvcfs", version_ch.collect())
	emit:
		version = version_ch
		region_vcf = genotyped_ch.region_vcf
}



process HAPCALLER_GVCF {
tag "bed:${file(row.bed).name} ref:${file(row.bam_reference).name} bam:${file(row.bam).name}"
cache 'lenient'
memory {task.attempt==1?'10G':'20G'}
errorStrategy 'retry'
maxRetries 5
cpus 1
afterScript 'rm -rf TMP'
input:
	val(meta)
	val(row)
output:
        tuple val(row),path("hapcaller.g.vcf.gz"),emit:bedvcf
        path("hapcaller.g.vcf.gz.tbi"),emit:index
        path("version.xml"),emit:version
script:
	def bam = row.bam?:""
	def sample = row.sample?:""
	def new_sample = row.new_sample?:""
	def extraHC = row.extraHC?:""
	def pedigree = row.pedigree?:""
	def mapq = row.mapq?:"-1"
	def bams = row.bams?:""
	def reference = row.reference?:""
	def bam_reference = row.bam_reference?:row.reference
	def bed = row.bed?:""
	def dbsnp = row.dbsnp?:""
"""
hostname 1>&2
${moduleLoad("gatk4 bcftools jvarkit")}

   set -o pipefail
   set -x
   mkdir TMP 

   if [ "${reference}" != "${bam_reference}" ] ; then
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/bedrenamechr.jar -R "${bam_reference}" --column 1 --convert SKIP "${bed}" > TMP/fixed.bed
	if [ ! -s TMP/fixed.bed ] ; then
		tail -n 1 "\${REF}.fai" | awk -F '\t' '{printf("%s\t0\t1\\n",\$1);}' > TMP/fixed.bed
   	fi 
   else
	ln -s "${bed}" TMP/fixed.bed
   fi


   gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" HaplotypeCaller \
     -I "${bam}" \
     -ERC GVCF \
     --seconds-between-progress-updates 600 \
     ${!isBlank(dbsnp) &&  reference.equals(bam_reference)?"--dbsnp ${dbsnp}":""}   \
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
	if [ ! -z "${dbsnp}" ] ; then
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
	val(meta)
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
	val(row)
output:
	path("bed_samplemap.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def contig = row.contig
	def blocksize= meta.blocksize?:"1mb"
	def mergesize = meta.mergesize?:"100"

if(parseBoolean(meta.use_whole_block))
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
	val(reference)
	path(pedigree)
	val(row)
output:
        tuple val("${row.interval}"),path("genotyped.bcf"),emit:region_vcf
        path("genotyped.bcf.csi"),emit:index
	path("version.xml"),emit:version
script:
	 def otherOpts1 =  gatkGetArgumentsForCombineGVCFs(meta.plus("reference":reference))
	 def otherOpts2 =  gatkGetArgumentsForGenotypeGVCF(meta.plus("reference":reference,"pedigree":pedigree))
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



