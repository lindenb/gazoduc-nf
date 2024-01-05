/*

Copyright (c) 2024 Pierre Lindenbaum

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
include { getVersionCmd;moduleLoad; getKeyValue; getModules; getBoolean; assertNotEmpty} from '../../modules/utils/functions.nf'
include { SAMTOOLS_SAMPLES as CASES_BAMS; SAMTOOLS_SAMPLES as CTRLS_BAMS} from '../../subworkflows/samtools/samtools.samples.03.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include { SCATTER_TO_BED } from '../../subworkflows/picard/picard.scatter2bed.nf'
include {SQRT_FILE} from '../../modules/utils/sqrt.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {CONCAT_FILES_01} from '../../modules/utils/concat.files.nf'
include {VCF_DUPHOLD_01} from '../duphold/vcf.duphold.01.nf'

workflow SMOOVE_SV_POPULATION_01 {
	take:
		meta
		genomeId
		cases
		controls
		exclude_bed
	main:

		version_ch = Channel.empty()

	
	
		each_cases = CASES_BAMS([:], cases).rows.map{T->T.plus("status":"case")}
		each_case_control_ch = each_cases
                //version_ch = version_ch.mix(each_case_control_ch.version)

                each_control = CTRLS_BAMS([:], controls).rows.map{T->T.plus("status":"control")}
		each_case_control_ch = each_case_control_ch.mix(each_control)



		
		each_sample_bam = each_case_control_ch.map{T->[T[0],T[1]]}
		
		if(exclude_bed.name.equals("NO_FILE")) {
			gaps_ch = SCATTER_TO_BED(["OUTPUT_TYPE":"N","MAX_TO_MERGE":"1"],params.genomes[genomeId].fasta )
			version_ch= version_ch.mix(gaps_ch.version)
			xbed = gaps_ch.bed
			}
		else
			{
			xbed = exclude_bed
			}

		img_ch = INSTALL_SMOOVE_IMAGE([:])
		version_ch= version_ch.mix(img_ch.version)

		call_ch = CALL_SMOOVE([:], genomeId, img_ch.smoove_img, xbed, each_cases)
		version_ch= version_ch.mix(call_ch.version.first())
	
		merge_ch = MERGE_ALL_SAMPLES([:], genomeId, img_ch.smoove_img, call_ch.vcf.collect())
		version_ch= version_ch.mix(merge_ch.version)

		gt_ch = GENOTYPE_BAM([:], genomeId, img_ch.smoove_img, merge_ch.vcf, each_sample_bam)
		version_ch= version_ch.mix(gt_ch.version.first())
		
		to_file_ch = COLLECT_TO_FILE_01([:], gt_ch.vcf.collect())
		sqrt_ch = SQRT_FILE([:], to_file_ch.output)
		version_ch= version_ch.mix(sqrt_ch.version)

		each_cluster = sqrt_ch.output.splitCsv(header: false,sep:',',strip:true).map{T->T[0]}
		
		paste01_ch = PASTE01([:], genomeId, img_ch.smoove_img, each_cluster)
		version_ch= version_ch.mix(paste01_ch.version)

		pasteall_ch = PASTE_ALL([:], genomeId, img_ch.smoove_img, paste01_ch.vcf.collect())
		version_ch= version_ch.mix(pasteall_ch.version)

		/** duphold is slow when many samples . */
		if(params.with_duphold) {
			all_bams_ch = CONCAT_FILES_01([:], Channel.from(cases,controls).collect())
			version_ch= version_ch.mix(all_bams_ch.version)

			duphold_ch = VCF_DUPHOLD_01([:], genomeId, all_bams_ch.output , pasteall_ch.vcf, file("NO_FILE"))
			version_ch= version_ch.mix(duphold_ch.version)
			
			final_vcf = duphold_ch.vcf
			final_idx = duphold_ch.index
			}
		else {
			final_vcf = pasteall_ch.vcf
			final_idx = pasteall_ch.index
		}

		version_ch = MERGE_VERSION("smoove", version_ch.collect())
	emit:
		version = version_ch.version
		vcf = final_vcf
		index = final_idx
	}


process INSTALL_SMOOVE_IMAGE {
	executor "local"
	input:
		val(meta)
	output:
		path("smoove.simg"),emit:smoove_img
		path("version.xml"),emit:version
	script:
		def img = params.smoove.singularity_image
	"""
	hostname 1>&2
	${moduleLoad("singularity/2.4.5")}
	set -x
	cp -v "${img}" "smoove.simg"

	singularity run smoove.simg smoove --help || true
	#######################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">test smoove img</entry>
		<entry key="smoove.image">${img}</entry>
		<entry key="smoove.version">\$(singularity run smoove.simg smoove --version)</entry>
	</properties>
	EOF
	"""
	}


/**
 *
 * environ 20-30 minutes pour un WGS cram
 *
 */
process CALL_SMOOVE {
tag "${row.sample}/${file(row.bam).name}"
cache "lenient"
memory "5g"
errorStrategy "finish"
afterScript  "rm -rf TMP TMP2"
cpus 1 /* can only parallelize up to 2 or 3 threads on a single-sample and it's most efficient to use 1 thread. */
input:
	val(meta)
	val(genomeId)
	path(img)
	val(exclude)
	val(row)
output:
	path("${row.sample}-smoove.genotyped.vcf.gz"),emit:vcf
	path("${row.sample}-smoove.genotyped.vcf.gz.tbi"),emit:tbi
	path("version.xml"),emit:version
script:
	def sample = row.sample
	def bam = row.bam
	def reference = params.genomes[genomeId].fasta
	def xbed = file(exclude)
	def ref = file(reference)
	def xbam = file(bam)
	if(xbam.getParent()==null) throw new IllegalArgumentException("no parent for ${bam}");
"""
	hostname 1>&2
	${moduleLoad("bcftools")}

	mkdir -p TMP
	# smoove will write to the system TMPDIR. For large cohorts, make sure to set this to something with a lot of space
	mkdir -p TMP2
	export TMPDIR=\${PWD}/TMP2

	singularity exec\
		--home \${PWD} \
		--bind ${xbed.getParent()}:/xdir \
		--bind ${xbam.getParent()}:/bamdir \
		--bind ${ref.getParent()}:/ref \
		--bind \${PWD}/TMP:/outdir \
		${img} \
		smoove call \
			--outdir /outdir \
			--exclude /xdir/${xbed.name} \
			--name ${sample} \
			--fasta /ref/${ref.name} \
			-p ${task.cpus} \
			--genotype \
			/bamdir/${xbam.name}

	mv TMP/${sample}-smoove.genotyped.vcf.gz ./
	bcftools index -t ${sample}-smoove.genotyped.vcf.gz
	
	rm -rf TMP

	#######################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">call with smoove</entry>
		<entry key="sample">${sample}</entry>
		<entry key="bam">${bam}</entry>
	</properties>
	EOF
"""
}

process MERGE_ALL_SAMPLES {
tag "N=${L.size()}"
cache "lenient"
memory "15g"
maxForks 1
cpus 16
afterScript  "rm -rf TMP TMP2"
input:
	val(meta)
	val(genomeId)
	val(img)
	val(L)
output:
	path("merged.sites.vcf.gz"),emit:vcf
	path("merged.sites.vcf.gz.tbi"),emit:tbi
	path("version.xml"),emit:version
when:
	L.size()>0

script:
	def reference = params.genomes[genomeId].fasta
	def ref = file(reference)
"""
	hostname 1>&2
	module load ${getModules("singularity bcftools")}

	mkdir -p TMP
	# smoove will write to the system TMPDIR. For large cohorts, make sure to set this to something with a lot of space
	mkdir -p TMP2
	export TMPDIR=\${PWD}/TMP2

	singularity exec\
		--home \${PWD} \
		${L.withIndex().collect{B,I -> " --bind "+ B.getParent()+":/data"+I}.join(" ")} \
		--bind ${ref.getParent()}:/ref \
		--bind \${PWD}/TMP:/outdir \
		${img} \
		smoove merge \
			--name merged \
			--outdir /outdir \
			--name merged \
			--fasta /ref/${ref.name} \
			${L.withIndex().collect{B,I -> "/data"+I+"/"+ B.name}.join(" ")}


	bcftools sort  --max-mem ${task.memory.giga}G  -T TMP -o merged.sites.vcf.gz -O z TMP/merged.sites.vcf.gz
	bcftools index -t -f merged.sites.vcf.gz


	#######################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">merge smoove calls</entry>
		<entry key="sample.count">${L.size()}</entry>
	</properties>
	EOF
"""
}


process GENOTYPE_BAM {
tag "${row.sample}/${file(row.bam).name}/${file(merged).name}"
cache "lenient"
errorStrategy "retry"
maxRetries 5
memory "10g"
afterScript  "rm -rf TMP TMP2"
cpus 1 /* can only parallelize up to 2 or 3 threads on a single-sample and it's most efficient to use 1 thread. */
input:
	val(meta)
	val(genomeId)
	val(img)
	val(merged)
	val(row)
output:
	path("${row.sample}-smoove.regenotyped.vcf.gz"),emit:vcf
	path("${row.sample}-smoove.regenotyped.vcf.gz.tbi"),emit:tbi
	path("version.xml"),emit:version
script:
	def reference = params.genomes[genomeId].fasta
	def sample = row.sample
	def bam = row.bam
	def ref = file(reference)
	def vcf0 = file(merged)
	def xbam = file(bam)
"""
	hostname 1>&2
	#module load singularity/2.4.5
	module load bcftools/0.0.0
	mkdir -p TMP
	# smoove will write to the system TMPDIR. For large cohorts, make sure to set this to something with a lot of space
	mkdir -p TMP2
	export TMPDIR=\${PWD}/TMP2

	singularity exec\
		--home \${PWD} \
		--bind ${xbam.getParent()}:/bamdir \
		--bind ${vcf0.getParent()}:/mergeddir \
		--bind ${ref.getParent()}:/ref \
		--bind \${PWD}/TMP:/outdir \
		${img} \
		smoove genotype -x  \
			--vcf /mergeddir/${vcf0.name} \
			--outdir /outdir \
			--name ${sample} \
			--fasta /ref/${ref.name} \
			-p ${task.cpus} \
			/bamdir/${xbam.name}

	bcftools sort  --max-mem ${task.memory.giga}G  -T TMP -o ${sample}-smoove.regenotyped.vcf.gz -O z TMP/${sample}-smoove.genotyped.vcf.gz
	bcftools index -t -f ${sample}-smoove.regenotyped.vcf.gz

	#######################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">genotype with smoove</entry>
		<entry key="sample">${sample}</entry>
		<entry key="bam">${bam}</entry>
		<entry key="merged.vcf">${merged}</entry>
	</properties>
	EOF

"""
}

process PASTE01 {
tag "${file(L).name}"
afterScript  "rm -rf TMP TMP2"
cache "lenient"
memory "10g"
input:
	val(val)
	val(genomeId)
	val(img)
	val(L)
output:
	path("paste.smoove.square.vcf.gz"),emit:vcf
	path("paste.smoove.square.vcf.gz.tbi"),emit:tbi
	path("version.xml"),emit:version
script:
	def reference = params.genomes[genomeId].fasta
"""
	hostname 1>&2

	module load bcftools/0.0.0
	#module load singularity/2.4.5
	mkdir -p TMP
	# smoove will write to the system TMPDIR. For large cohorts, make sure to set this to something with a lot of space
	mkdir -p TMP2
	export TMPDIR=\${PWD}/TMP2


	singularity exec\
		--home \${PWD} \
		--bind \${PWD}/TMP:/outdir \
		`xargs -a ${L} -L 1 dirname | awk '{printf(" --bind %s:/d%d ",\$0,NR);}'  ` \
		${img} \
		smoove paste --name paste \
			 `awk -F '/' '{printf(" /d%d/%s ",NR,\$NF);}' ${L}`


	bcftools index -t -f paste.smoove.square.vcf.gz

	#######################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">smoove paste</entry>
		<entry key="count">${L.size()}</entry>
	</properties>
	EOF

"""
}


process PASTE_ALL {
tag "N=${L.size()}"
cache "lenient"
afterScript  "rm -rf TMP TMP2"
memory "10g"
input:
	val(meta)
	val(reference)
	val(img)
	val(L)
output:
	path("${params.prefix?:""}smoove.bcf"),emit:vcf
	path("${params.prefix?:""}smoove.bcf.csi"),emit:index
	path("version.xml"),emit:version
script:
	def prefix = params.prefix?:""
	def gff = ""//gff3
	log.warn("JE SUPPRIME GFF POUR LE MOMENT. CA BUG POUR SOLENA OCT 2022")
"""
	hostname 1>&2
	${moduleLoad("bcftools picard")}
	set -x
	mkdir -p TMP TMP2
	# smoove will write to the system TMPDIR. For large cohorts, make sure to set this to something with a lot of space
	export TMPDIR=\${PWD}/TMP2

if [ "${L.size()}" -ne "1" ] ; then

cat << EOF > TMP/jeter.list
${L.join("\n")}
EOF

	#module load singularity/2.4.5
	singularity exec\
		--home \${PWD} \
		--bind \${PWD}/TMP:/outdir \
			`xargs -a TMP/jeter.list -L 1 dirname | awk '{printf(" --bind %s:/d%d ",\$0,NR);}'  ` \
		${img} \
		smoove paste --name paste  \
			 `awk -F '/' '{printf(" /d%d/%s ",NR,\$NF);}' TMP/jeter.list`

	mv paste.smoove.square.vcf.gz TMP/jeter1.vcf.gz

else

	bcftools view -O z -o TMP/jeter1.vcf.gz "${L[0]}"

fi


	if [ ! -z "${gff}" ] ; then
	
		singularity exec\
			--home \${PWD} \
			--bind \${PWD}/TMP:/outdir \
			${gff.isEmpty()?"":" --bind "+ file(gff).getParent()+":/gdir"} \
			--bind \${PWD}/TMP:/data1 \
			${img} \
			smoove annotate ${gff.isEmpty()?"":"--gff /gdir/"+file(gff).name} /data1/jeter1.vcf.gz > TMP/jeter2.vcf

		bcftools view -O z -o TMP/jeter2.vcf.gz TMP/jeter2.vcf

		mv TMP/jeter2.vcf.gz TMP/jeter1.vcf.gz
	fi


	# update dict
	java -jar -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP \${PICARD_JAR} UpdateVcfSequenceDictionary \
		I=TMP/jeter1.vcf.gz O=TMP/jeter2.vcf.gz SD=${reference}
	mv TMP/jeter2.vcf.gz TMP/jeter1.vcf.gz

	# sort vcf
	bcftools sort  --max-mem ${task.memory.giga}G  -T TMP -o "${prefix}smoove.bcf" -O b  TMP/jeter1.vcf.gz
	bcftools index -f "${prefix}smoove.bcf"

	#######################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">smoove paste all and annotate</entry>
		<entry key="gff">${gff}</entry>
		<entry key="count">${L.size()}</entry>
		<entry key="versions">${getVersionCmd("bcftools picard/UpdateVcfSequenceDictionary")}</entry>
	</properties>
	EOF
"""
}
