include {runOnComplete;moduleLoad} from '../../modules/utils/functions.nf'
include {SAMTOOLS_SAMPLES} from '../../subworkflows/samtools/samtools.samples.03.nf'


workflow  {
	CNVKIT([:], params.genomeId, file(params.bams), file(params.capture))
	}

workflow CNVKIT {
	take:
		meta
		genomeId
		bams
		capture
	main:
		sn_bam_ch = SAMTOOLS_SAMPLES([:],bams)

		refgene_ch = REFGENE([:], genomeId)
		access_ch = ACCESS([:], genomeId)
		target_ch = CNVKIT_TARGET([:], capture, refgene_ch.output)
		antitarget_ch = CNVKIT_ANTITARGET([:], capture, access_ch.output)
		refcnn_ch = CNVKIT_REFERENCE_NO_CONTROL([:], genomeId, target_ch.output, antitarget_ch.output)



		targets_ch = target_ch.output.map{T->["target",T]}.
			mix(antitarget_ch.output.map{T->["antitarget",T]}).
			map{T->[targettype:T[0],bed:T[1]]}
		



		rows1_ch = sn_bam_ch.rows.
			filter{T->T.genomeId.equals(genomeId)}.
			combine(targets_ch).
			map{T->T[0].plus(T[1])}




		coverage_ch = CNVKIT_COVERAGE([:], rows1_ch)
	
		covtarget_ch = coverage_ch.output.filter{T->T[0].targettype.equals("target")}.map{T->[T[0].sample,T[1]]}
		covtantiarget_ch = coverage_ch.output.filter{T->T[0].targettype.equals("antitarget")}.map{T->[T[0].sample,T[1]]}

		sn_target_antitarget_ch = covtarget_ch.join(covtantiarget_ch)

		/**
		OK
		here I haven't defined case or controls. I have two strategies
		 * I test the sample vs a dummy ref_cnn
		 * I test the sample vs a a ref_cnn of all other samples

		TODO add case/control list
		*/

		/** one sample vs the refcnn_ch.output */
		fix_single_ch = CNVKIT_FIX_NO_CONTROL([:],  refcnn_ch.output, sn_target_antitarget_ch )

		sn_vs_other_ch = sn_target_antitarget_ch.combine(sn_target_antitarget_ch).
			filter{T->!T[0].equals(T[3])}.
			map{T->[
				[sample:T[0],target:T[1],antitarget:T[2]] ,
				[T[4],T[5]] ]
				}.
			groupTuple().
			map{T->[T[0],T[1].flatten().sort()]}
	
		ref_vs_others_ch = CNVKIT_REFERENCE_SAMPLE_VS_OTHERS([:], genomeId, sn_vs_other_ch)

		fix_vs_others_ch  = CNVKIT_REFERENCE_FIX_VS_OTHERS([:], ref_vs_others_ch.output.map{T->T[0].plus(ref_cnn:T[1])} )


		MULTI_INTERSECT([:],
			fix_vs_others_ch.output.map{T->["vs_other",T[0].sample+","+T[4]]}.
			mix(fix_single_ch.output.map{T->["single",T[0]+","+T[4]]}).
			groupTuple()
			)
		
	}

process REFGENE {
	tag "${genomeId}"
	input:
		val(meta)
		val(genomeId)
	output:
		path("refFlat.${genomeId}.txt"),emit:output
	script:
		def genome = params.genomes[genomeId]
	"""
	${moduleLoad("jvarkit")}
	set -o pipefail
	wget -O - "http://hgdownload.soe.ucsc.edu/goldenPath/${genome.ucsc_name}/database/refFlat.txt.gz" |\\
		gunzip -c |\\
		java -jar  \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${genome.fasta}" --column 3 --convert SKIP > refFlat.${genomeId}.txt

	test -s "refFlat.${genomeId}.txt"
	"""
	}

process ACCESS {
	tag "${genomeId}"
        conda "${moduleDir}/../../conda/cnvkit.yml"
	input:
		val(meta)
		val(genomeId)
	output:
		path("access.${genomeId}.bed"),emit:output
	script:
		def genome = params.genomes[genomeId]
		def reference = genome.fasta
		def url = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/${genome.ucsc_name}-blacklist.v2.bed.gz"
	"""
	hostname 1>&2
	${moduleLoad("jvarkit")}
	set -o pipefail

	wget -O - "${url}" |\\
		gunzip -c |\\
		java -jar  \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${reference}" --convert SKIP > excludes.bed

	cnvkit.py access "${reference}" -x excludes.bed -o access.${genomeId}.bed
	"""
	}

process CNVKIT_TARGET {
        tag "${bed.name} ${refFlat.name}"
        conda "${moduleDir}/../../conda/cnvkit.yml"
        afterScript "rm -rf TMP"
	input:
		val(meta)
		path(bed)
		path(refFlat)
        output:
               	path("targets.bed"),emit:output
        script:
		def args = "--short-names --split --avg-size 200"		
        """
	hostname 1>&2
	mkdir -p TMP
        export TMPDIR=\${PWD}/TMP

        cnvkit.py target "${bed}" \\
                --annotate "${refFlat}" \\
		${args} \\
                -o targets.bed

        """
	}

/**
Given a target BED file that lists the chromosomal coordinates of the tiled regions used for targeted resequencing
 derive a BED file off-target/antitarget
**/
process CNVKIT_ANTITARGET {
        tag "${file(bed).name}"
        conda "${moduleDir}/../../conda/cnvkit.yml"
        afterScript "rm -rf TMP"
	input:
		val(meta)
		path(bed)
		path(access)
        output:
                path("antitargets.bed"),emit:output
        script:
		def args=""
        """
        mkdir TMP
        export TMPDIR=\${PWD}/TMP

        cnvkit.py antitarget \\
		${args} \\
		-g "${access}" \\
                -o antitargets.bed \\
                "${bed}"
        """
	}




/* ... for each samples */
process CNVKIT_COVERAGE {
	tag "${row.sample}"
	afterScript "rm -rf TMP"
        conda "${moduleDir}/../../conda/cnvkit.yml"
	input:
		val(meta)
		val(row)
	output:
		tuple val(row),path("${row.sample}.${row.targettype}coverage.cnn"),emit:output
	script:
		def genome = params.genomes[row.genomeId]
		def reference = genome.fasta
		def mapq = params.mapq?:30
	"""
	hostname 1>&2
	mkdir -p TMP
	export TMPDIR=\${PWD}/TMP

	cnvkit.py coverage \\
		-p ${task.cpus} \\
		-q ${mapq} \\
		--fasta "${reference}" \\
		-o ${row.sample}.${row.targettype}coverage.cnn \\
		'${row.bam}' '${row.bed}'
	"""
	}



process CNVKIT_REFERENCE_NO_CONTROL {
        conda "${moduleDir}/../../conda/cnvkit.yml"
	afterScript "rm -rf TMP"
	input:
		val(meta)
		val(genomeId)
		path(target)
		path(antitarget)
	output:
		path("${genomeId}.reference.cnn"),emit:output
	script:
		def genome = params.genomes[genomeId]
		def reference = genome.fasta
	"""
	mkdir TMP
	export TMPDIR=\${PWD}/TMP

   	 cnvkit.py reference \\
		-f "${reference}" \\
		-o "${genomeId}.reference.cnn" \\
		-t ${target} \\
		-a ${antitarget}

	"""
	}



/** For each tumor sample... */
process CNVKIT_FIX_NO_CONTROL {
	tag "${sample}"
        conda "${moduleDir}/../../conda/cnvkit.yml"
	afterScript "rm -rf TMP"
	input:
		val(meta)
		path(ref_cnn)
		tuple val(sample),path(target),path(antitarget)
	output:
		tuple 	val(sample),
			path("${sample}-scatter.pdf"),
			path("${sample}-diagram.pdf"),
			path("${sample}.cns.bed.gz"),
			path("${sample}.ncopy.bed.gz"),
			emit:output
	script:
	"""
	hostname >&2
	mkdir -p TMP
	export TMPDIR=\${PWD}/TMP

	## 23 juillet 2021 certains noeuds n'ont pas X11...
	unset DISPLAY

	cnvkit.py fix \\
		"${target}" \\
		"${antitarget}" \\
		"${ref_cnn}" \\
		-o ${sample}.cnr
	cnvkit.py segment ${sample}.cnr -o ${sample}.cns
	
	cnvkit.py scatter ${sample}.cnr -s ${sample}.cns -o ${sample}-scatter.pdf
	cnvkit.py diagram ${sample}.cnr -s ${sample}.cns -o ${sample}-diagram.pdf

	# Show estimated integer copy number of all regions
	cnvkit.py export bed "${sample}.cns" -o TMP/jeter.bed
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n TMP/jeter.bed > "${sample}.ncopy.bed"
	gzip --best "${sample}.ncopy.bed"
	
	sed 's/^chromosome/#chromosome/' "${sample}.cns" |\
		LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\
		gzip --best > "${sample}.cns.bed.gz"

	"""
	}


process CNVKIT_REFERENCE_SAMPLE_VS_OTHERS {
	tag "${key.sample} N=${L.size()}"
        conda "${moduleDir}/../../conda/cnvkit.yml"
	afterScript "rm -rf TMP"
	input:
		val(meta)
		val(genomeId)
		tuple val(key),val(L)
	output:
		tuple val(key),path("${key.sample}.reference.cnn"),emit:output
	script:
		def genome = params.genomes[genomeId]
		def reference = genome.fasta
	"""
	mkdir TMP
	export TMPDIR=\${PWD}/TMP

   	 cnvkit.py reference \\
		-f "${reference}" \\
		-o "${key.sample}.reference.cnn" \\
		${L.join(" ")}
	"""
	}

/** For each tumor sample... */
process CNVKIT_REFERENCE_FIX_VS_OTHERS {
	tag "${row.sample}"
        conda "${moduleDir}/../../conda/cnvkit.yml"
	afterScript "rm -rf TMP"
	input:
		val(meta)
		val(row)
	output:
		tuple 	val(row),
			path("${row.sample}-scatter.pdf"),
			path("${row.sample}-diagram.pdf"),
			path("${row.sample}.cns.bed.gz"),
			path("${row.sample}.ncopy.bed.gz"),
			emit:output
	script:
		def sample = row.sample
	"""
	hostname >&2
	mkdir -p TMP
	export TMPDIR=\${PWD}/TMP

	## 23 juillet 2021 certains noeuds n'ont pas X11...
	unset DISPLAY

	cnvkit.py fix \\
		"${row.target}" \\
		"${row.antitarget}" \\
		"${row.ref_cnn}" \\
		-o ${sample}.cnr
	cnvkit.py segment ${sample}.cnr -o ${sample}.cns
	
	cnvkit.py scatter ${sample}.cnr -s ${sample}.cns -o ${sample}-scatter.pdf
	cnvkit.py diagram ${sample}.cnr -s ${sample}.cns -o ${sample}-diagram.pdf

	# Show estimated integer copy number of all regions
	cnvkit.py export bed "${sample}.cns" -o TMP/jeter.bed
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n TMP/jeter.bed > "${sample}.ncopy.bed"
	gzip --best "${sample}.ncopy.bed"
	
	sed 's/^chromosome/#chromosome/' "${sample}.cns" |\
		LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\
		gzip --best > "${sample}.cns.bed.gz"

	"""
	}



process MULTI_INTERSECT {
tag "${type} N=${L.size()}"
input:
	val(meta)
	tuple val(type),val(L)
output:
	path("${params.prefix?:""}${type}.multiinter.bed.gz"),emit:output
script:
"""
${moduleLoad("bedtools")}

cat << EOF > jeter.csv
${L.join("\n")}
EOF

bedtools multiinter -header -names `cut -d, -f 1 jeter.csv` -i `cut -d, -f2 jeter.csv` > "${params.prefix?:""}${type}.multiinter.bed"
gzip --best "${params.prefix?:""}${type}.multiinter.bed"
rm jeter.csv
"""
}
