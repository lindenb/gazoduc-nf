/*

Copyright (c) 2022 Pierre Lindenbaum

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
include { getKeyValue; getModules; getBoolean} from '../../modules/utils/functions.nf'
include {DELLY2_RESOURCES} from './delly2.resources.nf' 
include {SAMTOOLS_CASES_CONTROLS_01} from '../samtools/samtools.cases.controls.01.nf'


workflow DELLY2_SV {
	take:
		meta
		reference
		cases
		controls
	main:
		version_ch = Channel.empty()
		rsrcr_ch = DELLY2_RESOURCES([:],reference)
		version_ch= version_ch.mix(rsrcr_ch.version)

		cases_controls_ch = SAMTOOLS_CASES_CONTROLS_01([:],reference,cases,controls)
		version_ch= version_ch.mix(cases_controls_ch.version)

		each_case_control_ch = cases_controls_ch.output.splitCsv(header:false,sep:'\t')
	
		each_cases = each_case_control_ch.filter{T->[2].equals("case")}.map{T->[T[0],T[1]]}

		delly_bcf = CALL_DELLY(meta, reference, rsrcr_ch.executable, rsrcr_ch.exclude, each_cases)
		version_ch= version_ch.mix(delly_bcf.version.first())

		merge_delly = 	MERGE_DELLY(meta, reference, rsrcr_ch.executable, delly_bcf.output.collect())
		version_ch= version_ch.mix(merge_delly.version)

		genotype_ch = GENOTYPE_DELLY(meta,reference, rsrcr_ch.executable, merge_delly.output, rsrcr_ch.exclude, each_case_control_ch.map{T->[T[0],T[1]]})
		version_ch= version_ch.mix(genotype_ch.version.first())

	
		merge_gt = MERGE_GENOTYPES(meta, genotype_ch.output.collect())
		version_ch= version_ch.mix(merge_gt.version)

		filter_delly = FILTER_DELLY(meta, rsrcr_ch.executable, cases_controls_ch.output, merge_gt.output ) 
		version_ch= version_ch.mix(filter_delly.version)
	emit:
		version = version_ch
		sv_vcf = filter_delly.output
	}


process CALL_DELLY {
    tag "${name}"
    memory '10 g'
    afterScript "rm -rf TMP"
    input:
	val(meta)
	val(reference)
	val(delly)
	val(exclude)
	tuple val(name),val(bam)
    output:
    	tuple val(name),path(bam),path("${name}.bcf"),emit:output
	path("version.xml"),emit:version
    script:
	"""
	hostname 1>&2
	mkdir -p TMP
	export TMPDIR=\${PWD}/TMP

	${delly} call --exclude "${exclude}" \
		--outfile "${name}.bcf" \
		--genome "${reference}" \
		"${bam}" 1>&2

	rm -rf TMP

	#######################################################################
	cat << EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">call delly</entry>
		<entry key="sample">${name}</entry>
		<entry key="bam">${bam}</entry>
	</properties>
	EOF
	"""
    }

// TODO merge many files: https://github.com/dellytools/delly/issues/158
process MERGE_DELLY {
    tag "N=${bcfs.size()}"
    memory '10 g'
    input:
	val(meta)
	val(reference)
	path(delly)
	val(bcfs)
    output:
	path("merged.bcf"),emit:output
	path("version.xml"),emit:version
    script:
	def bnd = getBoolean(meta,"bnd")
    """
    hostname 1>&2
    module load getModules("jvarkit bcftools")}
    export LC_ALL=C
    mkdir TMP
    export TMPDIR=\${PWD}/TMP

# see https://github.com/dellytools/delly/issues/158
cat << EOF > TMP/jeter.tsv
${bcfs.join("\n")}
EOF
    
    ${delly} merge -o merged.bcf TMP/jeter.tsv 1>&2

	if [ "${bnd?"Y":"N"}" == "N" ] ; then
                bcftools view -e 'INFO/SVTYPE="BND"' -O b -o  jeter.bcf merged.bcf
		mv jeter.bcf merged.bcf
		bcftools index -f merged.bcf
	fi


	#######################################################################
	cat << EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">merge delly</entry>
		<entry key="keep bnd">${bnd}</entry>
	</properties>
	EOF

    """
    }

process GENOTYPE_DELLY {
    tag "${name}"
    cache 'lenient'
    errorStrategy 'finish'
    afterScript 'rm -rf TMP'
    memory '10 g'
    input:
	val(meta)
	val(reference)
	val(delly)
        val(merged)
	val(exclude)
        tuple val(name),val(bam)
    output:
        path("genotyped.${name}.bcf"),emit:output
	path("version.xml"),emit:version
    script:
    """
    hostname 1>&2
    mkdir -p TMP
    export TMPDIR=\${PWD}/TMP

    ${delly} call --vcffile "${merged}" \
		--exclude "${exclude}" \
		--outfile "TMP/jeter.bcf" \
		--genome "${reference}" \
		${bam} 1>&2

    mv -v TMP/jeter.bcf "genotyped.${name}.bcf"
    mv -v TMP/jeter.bcf.csi "genotyped.${name}.bcf.csi"

	#######################################################################
	cat << EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">genotype delly</entry>
		<entry key="sample">${name}</entry>
		<entry key="bam">${bam}</entry>
		<entry key="merged">${merged}</entry>
	</properties>
	EOF
    """
    }


process MERGE_GENOTYPES {
    tag "N=${bcfs.size()}"
    cache 'lenient'
    memory "20g"
    input:
	val(meta)
	val(bcfs)
    output:
	path("merged.gt.bcf"),emit:output
	path("version.xml"),emit:version
	path("merged.gt.bcf.csi")
    script:
    """
    hostname 1>&2
    module load ${getModules("bcftools")}

cat << EOF > jeter.list
${bcfs.join("\n")}
EOF

    bcftools merge -m id -O b -o merged.gt.bcf --file-list jeter.list
    bcftools index --csi merged.gt.bcf 
    rm jeter.list


	#######################################################################
	cat << EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">merge genotype delly</entry>
		<entry key="number of files">${L.size()}</entry>
	</properties>
	EOF

    """
    }

process FILTER_DELLY {
    tag "filter"
    cache 'lenient'
    memory "5g"
    input:
        val(meta)
        val(delly)
	val(cases_ctrl_list)
	val(merged)
    output:
	path("${prefix}bcf"),emit:output
	path("version.xml"),emit:version
	path("${prefix}bcf.csi")
    script:
	prefix = getKeyValue(meta,"prefix","")
    """
    export LC_ALL=C
    module load ${getModules("bcftools")}
    ${delly} filter -f germline  -o jeter.bcf "${merged}" 1>&2

    bcftools sort --max-mem "${task.memory.giga}G" -T . -O v -o "jeter1.vcf" jeter.bcf


    awk -F '\t' '(\$3=="case") {printf("%s\\n",\$1);}' "${cases_ctrl_list}" | sort | uniq > jeter.cases.txt
    awk -F '\t' '(\$3=="control") {printf("%s\\n",\$1);}' "${cases_ctrl_list}"| sort | uniq > jeter.ctrls.txt

    if [ ! -s "jeter.cases.txt" ] && [ ! -s "jeter.ctrls.txt"	] ; then
	# rajoute mais pas teste
	bcftools +contrast \
		-0 jeter.ctrls.txt \
		-1 jeter.cases.txt \
		-a PASSOC,FASSOC,NASSOC,NOVELAL,NOVELGT -O v -o jeter2.vcf jeter1.vcf
	mv jeter2.vcf jeter1.vcf
    fi
    rm jeter.cases.txt jeter.ctrls.txt
    

    bcftools view -O b -o "${prefix}bcf" jeter1.vcf
    bcftools index "${prefix}bcf"

    rm jeter.bcf jeter1.vcf

	#######################################################################
	cat << EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">delly filter</entry>
	</properties>
	EOF

    """
    }

/********************

process callCNV {
    tag "${name}"
    cache 'lenient'
    errorStrategy 'finish'
    memory '10 g'
    input:
	val(mappability) from mappability_out
	path delly from executable
	tuple name,bam,sv from delly_bcf3
    output:
    	tuple name,bam,path("${name}.cnv.bcf") into (delly_cnv1,delly_cnv2)
    when:
	!mappability.isEmpty() && params.cnv==true
    script:
	"""
	hostname 1>&2
	mkdir -p TMP
	export TMPDIR=\${PWD}/TMP
	${delly} cnv \
		--outfile "${name}.cnv.bcf" \
		--mappability "${mappability}" \
		--genome "${params.reference}" \
		--svfile "${sv.toRealPath()}" \
		"${bam}" 1>&2
	rm -rf TMP
	"""
    }

process mergeCnv {
    tag "N=${bcfs.size()}"
    cache 'lenient'
    memory '20 g'
    input:
	path delly from executable
	val bcfs from delly_cnv1.map{T->T[2]}.collect()
    output:
	path("merged.bcf") into merged_cnv
    when:
	params.cnv==true
    script:
    """
    hostname 1>&2
    export LC_ALL=C

cat << EOF > jeter.tsv
${bcfs.join("\n")}
EOF

    ${delly} merge -e -p -o merged.bcf  -m 1000 -n 10000000 jeter.tsv 1>&2

rm jeter.tsv

        if [ -f "${params.commonBed}" ] ; then
		bcftools view -O v merged.bcf |\
                java -Xmx${task.memory.giga}G -jar ${jvarkit("vcfbed")} \
                        --bed '${params.commonBed}'  \
                        --fast \
                        --min-overlap-bed-fraction ${params.fraction} \
                        --min-overlap-vcf-fraction ${params.fraction} \
                        -T COMMON  > jeter.vcf

                bcftools view -e 'INFO/COMMON!=""' -O b -o  merged.bcf jeter.vcf
		bcftools index -f merged.bcf
		rm jeter.vcf
        fi

    """
    }


process genotypeCNV {
    tag "${name}"
    cache 'lenient'
    errorStrategy 'finish'
    afterScript 'rm -f jeter.bcf jeter.bcf.csi'
    memory '10 g'
    input:
	val(mappability) from mappability_out
	path delly from executable
        path merged from merged_cnv
        tuple name,bam from sample_bam4
    output:
        path("genotyped.${name}.cnv.bcf") into gt_cnv_bcf
    when:
	params.cnv==true
    script:
    """
    hostname 1>&2
    mkdir -p TMP
    export TMPDIR=\${PWD}/TMP
    # je baisse min clip pour que le genotypage soit plus sensible
    ${delly} cnv  \
		--segmentation \
                --vcffile "${merged}" \
		--outfile "jeter.bcf" \
		--mappability "${mappability}" \
		--genome "${params.reference}" \
		${bam} 1>&2

    mv -v jeter.bcf "genotyped.${name}.cnv.bcf"
    mv -v jeter.bcf.csi "genotyped.${name}.cnv.bcf.csi"
    rm -rf TMP
    """
    }

process mergeCNV {
    tag "merge GT N=${bcfs.size()}"
    cache 'lenient'
    input:
	val bcfs from gt_cnv_bcf.collect()
    output:
	path("merged.gt.bcf") into mergedgt_cnv_bcf
	path("merged.gt.bcf.csi")
    when:
	params.cnv
    script:
    """
    module load bcftools/0.0.0
    bcftools merge -m id -O b -o merged.gt.bcf ${bcfs.join(" ")}
    bcftools index --csi merged.gt.bcf 
    """
    }

process classifyCNV {
    cache 'lenient'
    publishDir params.publishDir , mode: 'copy', overwrite: true
    afterScript "rm -f jeter.bcf jeter1.vcf jeter2.vcf"
    cpus 5
    memory "5g"
    input:
	path delly from executable
	path merged from mergedgt_cnv_bcf
    output:
	path("${params.prefix}cnv.bcf") into delly_cnv_output
	path("${params.prefix}cnv.bcf.csi")
    when:
	params.cnv
    script:
    """
    export LC_ALL=C
    module load bcftools/0.0.0
    ${delly} classify -f germline  -o jeter.bcf "${merged.toRealPath()}" 1>&2

    bcftools sort --max-mem ${task.memory.giga}G -T . -O v -o "jeter1.vcf" jeter.bcf

    if [ "${params.gtf}" != "" ] ; then
	    java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar ${jvarkit("svpredictions")} --max-genes 30  --gtf "${params.gtf}" jeter1.vcf >  jeter2.vcf
	    mv jeter2.vcf jeter1.vcf
    fi

    bcftools view -O b -o "${params.prefix}cnv.bcf" jeter1.vcf
    bcftools index "${params.prefix}cnv.bcf"

    """
    }

*/
