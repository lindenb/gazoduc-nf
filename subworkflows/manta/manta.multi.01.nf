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
include {moduleLoad;getVersionCmd} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'

workflow MANTA_MULTI_SV01 {
	take:
		meta
		genomeId
		bams
		bed
	main:
		version_ch = Channel.empty()

		config_ch = MANTA_MULTI([:], genomeId , bams,bed)
		version_ch= version_ch.mix(config_ch.version)
	
		
		version_ch = MERGE_VERSION("manta",version_ch.collect())
		html = VERSION_TO_HTML(version_ch.version)
	}

process MANTA_MULTI {
    tag "${bams.name}"
    afterScript "rm -rf TMP/workspace"
    cpus 16
    memory "20g"
    input:
	val(meta)
	val(genomeId)
	path(bams)
	path(bed)

    output:
	path("*.bcf"),emit:output
	path("*.csi")
	path("version.xml"),emit:version

    script:
	def reference = params.genomes[genomeId].fasta
	"""
	hostname 1>&2
	${moduleLoad("manta htslib bcftools samtools")}
	mkdir -p TMP

	if ${bed.name.equals("NO_FILE")} ; then

		awk -F '\t' '(\$1 ~ /^(chr)?[0-9XY]+\$/ ) {printf("%s\t0\t%s\\n",\$1,\$2)}' '${reference}.fai' |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n > manta.bed
	else
		${bed.name.endsWith(".gz")?"gunzip -c":"cat"} ${bed} |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n > manta.bed
	
	fi

	bgzip manta.bed
	tabix -fp bed manta.bed.gz

	configManta.py  `awk '{printf("--bam %s ",\$0);}' "${bams}" ` \
		--referenceFasta "${reference}" \
		--callRegions manta.bed.gz \
		--runDir "TMP"


	./TMP/runWorkflow.py --quiet -m local -j ${task.cpus}
	rm -rf TMP/workspace

        # convert BND TO INVERSIONS (added 20230115 but not tested)
        DIPLOID=`find ./TMP -type f -name "diploidSV.vcf.gz"`
        test ! -z "\${DIPLOID}"
        \$(ls \$( dirname \$(which configManta.py) )/../share/manta*/libexec/convertInversion.py)  `which samtools` "${reference}" "\${DIPLOID}" | bcftools sort -T TMP -O z -o TMP/jeter.vcf.gz

        bcftools index -t TMP/jeter.vcf.gz

        mv -v TMP/jeter.vcf.gz "\${DIPLOID}"
        mv -v TMP/jeter.vcf.gz.tbi "\${DIPLOID}.tbi"

	
	for X in diploidSV candidateSmallIndels candidateSV
	do
		bcftools view -O b -o "${params.prefix?:""}\${X}.bcf" "\${PWD}/TMP/results/variants/\${X}.vcf.gz"
		bcftools index "${params.prefix?:""}\${X}.bcf"
	done


#################################################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="versions">${getVersionCmd("bcftools")}</entry>
	<entry key="manta.version">\$(configManta.py --version)</entry>
</properties>
EOF
	"""
	}


