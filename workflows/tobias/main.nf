
workflow {
	TOBIAS(params.fasta, file(params.blacklist),file(params.samplesheet))
	}

def TOBIAS_CONTAINER = "https://depot.galaxyproject.org/singularity/tobias:0.13.3--py37h37892f8_0"

workflow TOBIAS {
	take:
		fasta
		blacklist
		samplesheet
	main:
		ref_ch = Channel.of([file(fasta),file(fasta+".fai")])

		conditions_ch = Channel.of(params.conditions.split("[,; ]")).
			filter{!it.trim().isEmpty()}
		

		samples_ch = Channel.fromPath(samplesheet).splitCsv(header:true,sep:'\t')


		/** extract all distincts condition from samplesheet . At the end we have a channel of [CONDITION_NAME,CONDITION_VALUE ] */
		cond_value_ch = samples_ch.combine(conditions_ch).
			filter{ it[0].containsKey(it[1])  }.
			filter{ it[0].containsKey(it[1])  }.
			map{ [ it[1] , it[0].get(it[1]) ] }.
			filter{!it[1].trim().isEmpty()}.
			unique()


		

		/** extract all pairs of comparaitions. At the end we have a channel of [CONDITION_NAME,VALUE1,VALUE2] */
		cond_v1_vs_v2_ch = cond_value_ch.combine(cond_value_ch).
			filter{it[0].equals(it[2]) && it[1].compareTo(it[3])<0}.
			map{[it[0],it[1],it[3]]}

		ch1 = samples_ch.map{
			T->[
				T.id,
				file(T.bam),
				file(T.containsKey("bai")?T.bai:T.bam+".bai"),
				file(T.peaks)
				]
			}


		

		correct_ch = ATACORRECT(blacklist, ref_ch.combine(ch1))

		ch2 = correct_ch.combine(samples_ch).
			filter{it[0].equals(it[2].id)}.
			map{
			T->[
				T[0],//id
				T[1],//corrected
				file(T[2].peaks)
				]
			}

		footprint_ch = FOOTPRINT_SCORE(ch2)
		
		ch3 = footprint_ch.combine(samples_ch).
			filter{it[0].equals(it[2].id)}


		
		jaspar1_ch = DOWNLOAD_JASPAR()
		jaspar2_ch = FORMAT_MOTIFS(jaspar1_ch.output)
		
		peak_header_ch = PEAK_HEADER()


		ch3 = footprint_ch.combine(samples_ch).
			filter{it[0].equals(it[2].id)}.
			map{
			T->[
				T[2],//sample
				T[1],//footprint
				file(T[2].peaks)
				]
			}

		ch4  = ch3.map{T->[T[0].id,T[1],T[2]]}

/*
		BIND_DETECT_SINGLE(
			jaspar2_ch.output,
			peak_header_ch.output,
			ref_ch.combine(ch4)
			)
*/
		ch5 = cond_value_ch.combine(ch3).
			filter{T->T[2].get(T[0]).equals(T[1])}.
			map{T->[
				condition_key : T[0],
				condition_value : T[1],
				id : T[2].id,
				footprint : T[3],
				peak : T[4]
				]
			}
		
		ch6 = ch5.combine(ch5).
			filter{T->T[0].condition_key.equals(T[1].condition_key)}.
			filter{T->T[0].condition_value.compareTo(T[1].condition_value)<0}.
			map{T->[
				T[0].condition_key, // condition

				T[0].condition_value, // cond1
				T[0].id,
				T[0].footprint,
				T[0].peak,
				
				T[1].condition_value, // cond1
				T[1].id,
				T[1].footprint,
				T[1].peak
				]}

		BIND_DETECT(jaspar2_ch.output,
                        peak_header_ch.output,
                        ref_ch.combine(ch6)
                        )

	}


/**

Similar to other enzymes used in chromatin accessibility assays (e.g. DNaseI for DNase-seq), the Tn5 transposase harbours an inherent insertion bias.
This means that while we assume insertion frequency to be driven by accessibility alone, the local cutting pattern is largely driven by the underlying sequence.
This interferes with footprinting analysis, and should therefore be corrected.

*/
process ATACORRECT {
afterScript "rm -rf TMP"
tag "${id}"
container  "${TOBIAS_CONTAINER}"
input:
	path(blacklist)
	tuple path(fasta),path(fai),val(id),path(bam),path(bai),path(peaks)
output:
	tuple val(id),path("${bam.getBaseName()}_*"),emit: output
script:
"""
hostname 1>&2
mkdir -p TMP

TOBIAS ATACorrect \\
	--bam "${bam}" \\
	--genome "${fasta}" \\
	--peaks "${peaks}" \\
	--blacklist "${blacklist}" \\
	--outdir TMP \\
	--cores "${task.cpus}" 2>&1 >  "${id}.atacorrect.log"

mv -v TMP/${bam.getBaseName()}_* ./
"""
}




/**

The main task in footprinting is to identify regions of protein binding across the genome.
Using single basepair cutsite tracks (as produced by ATACorrect),
TOBIAS ScoreBigwig is used to calculate a continuous footprinting score across regions.

*/
process FOOTPRINT_SCORE {
afterScript "rm -rf TMP"
tag "${id}"
container  "${TOBIAS_CONTAINER}"
input:
	tuple val(id),path(corrected),path(peaks)
output:
	tuple val(id),path("${id}.footprint.bw"),emit:output
script:
"""
hostname 1>&2
mkdir -p TMP

TOBIAS FootprintScores \\
	--signal *_corrected.bw \\
	--regions "${peaks}" \\
	--output TMP/${id}.footprint.bw \\
	--score footprint \\
	--cores "${task.cpus}" 2>&1 >  "${id}.footprint.log"

mv -v TMP/${id}.footprint.bw ./
"""
}


process DOWNLOAD_JASPAR {
executor "local"
afterScript "rm -rf TMP"
output:
	path("jaspar"),emit:output
script:
"""
mkdir -p TMP
wget -O TMP/jeter.zip  "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_non-redundant_pfms_jaspar.zip"
(cd TMP && unzip -o jeter.zip && rm -vf jeter.zip)
mv -v TMP jaspar
"""
}

process FORMAT_MOTIFS {
container  "${TOBIAS_CONTAINER}"
input:
	path(jaspar)
output:
	path("joined_motifs.jaspar"),emit:output
script:
"""
mkdir -p TMP
TOBIAS FormatMotifs \\
	--input ${jaspar}/* \\
	--task join \\
	--output TMP/joined_motifs.jaspar

if ${!params.tfb_regex.isEmpty()} ; then

	cat "TMP/joined_motifs.jaspar" | paste -d '|' - - - - - | grep -E '${params.tfb_regex}' | tr '|' '\\n' > TMP/tmp.jaspar

	mv -v TMP/tmp.jaspar TMP/joined_motifs.jaspar

fi

test -s TMP/joined_motifs.jaspar

mv TMP/joined_motifs.jaspar ./
"""
}

process PEAK_HEADER {
executor "local"
output:
	path("peak.header"),emit:output
script:
if(params.peakType.equals("narrow"))
"""

cat << EOF > peak.header
chrom
start
end
name
score
strand
signalValue
pValue
qValue
peak
EOF

"""
else
"""
false
"""
}

/**

BINDetect can also be used to predict binding for a single condition (which will turn off any estimation of differential binding):

*/
process BIND_DETECT_SINGLE {
tag "${id}"
afterScript "rm -rf TMP"
container  "${TOBIAS_CONTAINER}"
input:
	path(jaspar)
	path(peak_header)
	tuple path(fasta),
		path(fai),
		val(id),
		path(footprint),
		path(peaks)
output:
	tuple val(id),path("${id}.tf.txt"),path("${id}.OUTPUT"),emit:output
script:
"""
hostname 1>&2
mkdir -p TMP


TOBIAS BINDetect \\
	--motifs "${jaspar}" \\
	--signals ${footprint} \\
	--genome ${fasta} \\
	--peaks ${peaks} \\
	--peak_header ${peak_header} \\
	--outdir TMP \\
	--cond_names ${id} \\
	--cores ${task.cpus} 2> ${id}.bindetect.log

mv  "TMP" "${id}.OUTPUT"
find "${id}.OUTPUT" -maxdepth 1 -type d | awk -F '/' '(\$2!="") {print \$2}' > "${id}.tf.txt"
"""
}



process BIND_DETECT {
tag "${condition} ${cond1}/${id1} vs ${cond2}/${id2}"
afterScript "rm -rf TMP"
container  "${TOBIAS_CONTAINER}"
cpus 5
input:
	path(jaspar)
	path(peak_header)
	tuple path(fasta),
		path(fai),
		val(condition),
		val(cond1),
		val(id1),
		path(footprint1),
		path(peak1),
		val(cond2),
		val(id2),
		path(footprint2),
		path(peak2)

output:
	tuple val(condition),val(cond1),val(id1),val(cond2),val(id2),path("${id1}.${cond1}.${id2}.${cond2}.tf.txt"),path("${id1}.${cond1}.${id2}.${cond2}.OUTPUT"),emit:output
script:
"""
hostname 1>&2
mkdir -p TMP

cat ${peak1} ${peak2} > TMP/merged.bed


TOBIAS BINDetect \\
	--motifs "${jaspar}" \\
	--signals ${footprint1} ${footprint2} \\
	--genome ${fasta} \\
	--peaks TMP/merged.bed \\
	--peak_header ${peak_header} \\
	--outdir TMP \\
	--cond_names ${condition}_${cond1} ${condition}_${cond2} \\
	--cores ${task.cpus} 2> ${id1}.${cond1}.${id2}.${cond2}.bindetect.log

mv  "TMP" "${id1}.${cond1}.${id2}.${cond2}.OUTPUT"
find "${id1}.${cond1}.${id2}.${cond2}.OUTPUT" -maxdepth 1 -type d | awk -F '/' '(\$2!="") {print \$2}' > "${id1}.${cond1}.${id2}.${cond2}.tf.txt"
"""
}

