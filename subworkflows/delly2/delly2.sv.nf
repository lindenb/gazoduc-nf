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
include {getVersionCmd;moduleLoad; getBoolean} from '../../modules/utils/functions.nf'
include {DELLY2_RESOURCES} from './delly2.resources.nf' 
include { SAMTOOLS_SAMPLES as CASES_BAMS; SAMTOOLS_SAMPLES as CTRLS_BAMS} from '../../subworkflows/samtools/samtools.samples.03.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'

workflow DELLY2_SV {
	take:
		meta
		genomeId
		cases
		controls
		genotype_vcf /** SV to genotype or NO_FILE */
		user_exclude_bed /** cutsom bed exclude or NO_FILE */
		cnv_bed /** bed for CNV or NO_FILE */
	main:
		version_ch = Channel.empty()
		rsrcr_ch = DELLY2_RESOURCES([:], genomeId)
		version_ch= version_ch.mix(rsrcr_ch.version)


		each_case_control_ch = CASES_BAMS([:], cases).rows.map{T->T.plus("status":"case")}
                //version_ch = version_ch.mix(each_case_control_ch.version)

                if(!controls.name.equals("NO_FILE")) {
                        ctrls_ch = CTRLS_BAMS([:], controls).rows.map{T->T.plus("status":"control")}
                        //version_ch = version_ch.mix(ctrls_ch.version)
			each_case_control_ch = each_case_control_ch.mix(ctrls_ch)
                        }

		each_case_control_ch.combine(each_case_control_ch).
			filter{T->T[0].sample.equals(T[1].sample) && !T[0].bam.equals(T[1].bam)}.
			map{T->{
				throw new IllegalStateException("duplicate sample ${T}");
				}}

		each_case_control_ch = each_case_control_ch.
			filter{T->T.genomeId.equals(genomeId)}.
			map{T->[T.sample,T.bam,T.status]}


		each_sample_bam = each_case_control_ch.map{T->[T[0],T[1]]}

		if(user_exclude_bed.name.equals("NO_FILE")) {
			exclude_bed = rsrcr_ch.exclude
			}
		else {
			exclude_bed  = user_exclude_bed
			}


		/* run SV discovery */
		if( genotype_vcf.name.equals("NO_FILE")) {
			each_cases = each_case_control_ch.filter{T->T[2].equals("case")}.map{T->[T[0],T[1]]}
		
		
			delly_bcf = CALL_DELLY([:], genomeId, rsrcr_ch.executable, exclude_bed , each_cases)
			version_ch= version_ch.mix(delly_bcf.version.first())

			merge_delly = 	MERGE_DELLY([:], genomeId, rsrcr_ch.executable, delly_bcf.output.map{T->T[2]}.collect())
			version_ch= version_ch.mix(merge_delly.version)
	
			call_vcf = merge_delly.output	
			}
		else
			{
			call_vcf = genotype_vcf
			}

		genotype_ch = GENOTYPE_DELLY([:],genomeId, rsrcr_ch.executable, call_vcf, exclude_bed, each_sample_bam)
		version_ch= version_ch.mix(genotype_ch.version.first())

	
		merge_gt = MERGE_GENOTYPES([:], genotype_ch.output.collect())
		version_ch= version_ch.mix(merge_gt.version)

		filter_delly = FILTER_DELLY([:], rsrcr_ch.executable, each_case_control_ch.collect() , merge_gt.output ) 
		version_ch= version_ch.mix(filter_delly.version)

		if(params.cnv==true) {
			if( genotype_vcf.name.equals("NO_FILE")) {

				cnv_bcf = CALL_CNV([:], genomeId, rsrcr_ch.executable,  rsrcr_ch.mappability, cnv_bed, each_cases.combine(filter_delly.output) )
				version_ch = version_ch.mix(cnv_bcf.version.first())

				cnv_merge = MERGE_CNV([:], rsrcr_ch.executable,cnv_bcf.output.map{T->T[2]}.collect())
				version_ch = version_ch.mix(cnv_merge.version)
				
				call_cnv = cnv_merge.output
				}
			else	
				{
				call_cnv = genotype_vcf
				}

			cnv_genotype = GENOTYPE_CNV([:], genomeId, rsrcr_ch.executable, call_cnv, rsrcr_ch.mappability, cnv_bed,  each_sample_bam)
			version_ch = version_ch.mix(cnv_genotype.version.first())
	
			cnv_gt_merge = MERGE_CNV_GENOTYPED([:], cnv_genotype.output.collect())			
			version_ch = version_ch.mix(cnv_gt_merge.version)

			classify = CLASSIFY_CNV([:], rsrcr_ch.executable, cnv_gt_merge.output)
			version_ch = version_ch.mix(classify.version)

			cnv_output = classify.output
			cnv_index = classify.index
			} 
		else	{
			cnv_output = Channel.empty()
			cnv_index = Channel.empty()
			}
		version_merge = MERGE_VERSION("delly2",version_ch.collect())
	emit:
		version = version_merge.version
		sv_vcf = filter_delly.output
		sv_vcf_index = filter_delly.index
		cnv_vcf = cnv_output
		cnv_vcf_index = cnv_index
	}


process CALL_DELLY {
    tag "${name}"
    memory "10g"
    afterScript "rm -rf TMP"
    input:
	val(meta)
	val(genomeId)
	path(delly)
	val(exclude)
	tuple val(name),val(bam)
    output:
    	tuple val(name),val(bam),path("${name}.bcf"),emit:output
	path("version.xml"),emit:version
    script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
	"""
	hostname 1>&2
	mkdir -p TMP
	export TMPDIR=\${PWD}/TMP
	export PATH=\${PWD}:\${PATH}

	delly call --exclude "${exclude}" \
		--outfile "TMP/${name}.bcf" \
		--genome "${reference}" \
		"${bam}" 1>&2

	mv "TMP/${name}.bcf" ./

	#######################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">call delly</entry>
		<entry key="sample">${name}</entry>
		<entry key="bam">${bam}</entry>
		<entry key="versions">${getVersionCmd("delly")}</entry>
	</properties>
	EOF
	"""
    }

//  merge many files: https://github.com/dellytools/delly/issues/158
process MERGE_DELLY {
    tag "N=${bcfs.size()}"
    afterScript "rm -rf TMP"
    memory "20g"
    input:
	val(meta)
	val(genomeId)
	path(delly)
	val(bcfs)
    output:
	path("merged.bcf"),emit:output
	path("version.xml"),emit:version
	path("merged.bcf.csi")
    script:
        def genome = params.genomes[genomeId]
        def reference =	genome.fasta
	def bnd = params.bnd
    """
    hostname 1>&2
    ${moduleLoad("bcftools")}
    export LC_ALL=C
    mkdir TMP
    export TMPDIR=\${PWD}/TMP
    export PATH=\${PWD}:\${PATH}

# see https://github.com/dellytools/delly/issues/158
cat << EOF > TMP/jeter.tsv
${bcfs.join("\n")}
EOF
    
    delly merge -o TMP/merged.bcf TMP/jeter.tsv 1>&2

	if ${!bnd} ; then
                bcftools view -e 'INFO/SVTYPE="BND"' -O b -o  merged.bcf TMP/merged.bcf
		bcftools index -f merged.bcf
	else
		mv TMP/merged.bcf ./
		mv TMP/merged.bcf.csi ./
	fi


	#######################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">merge delly</entry>
		<entry key="keep bnd">${bnd}</entry>
		<entry key="versions">${getVersionCmd("bcftools")}</entry>
	</properties>
	EOF

    """
    }

process GENOTYPE_DELLY {
    tag "${name}"
    cache 'lenient'
    errorStrategy 'finish'
    afterScript 'rm -rf TMP'
    memory "15g"
    input:
	val(meta)
	val(genomeId)
	path(delly)
        val(merged)
	val(exclude)
        tuple val(name),val(bam)
    output:
        path("genotyped.${name}.bcf"),emit:output
	path("version.xml"),emit:version
    script:
        def genome = params.genomes[genomeId]
        def reference =	genome.fasta
    """
    hostname 1>&2
    mkdir -p TMP
    export TMPDIR=\${PWD}/TMP
    export PATH=\${PWD}:\${PATH}

    delly call --vcffile "${merged}" \
		--exclude "${exclude}" \
		--outfile "TMP/jeter.bcf" \
		--genome "${reference}" \
		"${bam}" 1>&2

    mv -v TMP/jeter.bcf "genotyped.${name}.bcf"
    mv -v TMP/jeter.bcf.csi "genotyped.${name}.bcf.csi"

	#######################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">genotype delly</entry>
		<entry key="sample">${name}</entry>
		<entry key="bam">${bam}</entry>
		<entry key="merged">${merged}</entry>
		<entry key="versions">${getVersionCmd("delly")}</entry>
	</properties>
	EOF
    """
    }


process MERGE_GENOTYPES {
    tag "N=${bcfs.size()}"
    afterScript "rm -rf TMP"
    label "process_high"
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

   if(bcfs.size() < 1000 )
    """
    hostname 1>&2
    ${moduleLoad("bcftools")}
    ulimit -s unlimited

cat << EOF > jeter.list
${bcfs.join("\n")}
EOF

    bcftools merge -m id -O b -o merged.gt.bcf --file-list jeter.list
    bcftools index --csi merged.gt.bcf 
    rm jeter.list


	#######################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">merge genotype delly</entry>
		<entry key="versions">${getVersionCmd("bcftools")}</entry>
		<entry key="number of files">${bcfs.size()}</entry>
	</properties>
	EOF

    """
else
    """
    hostname 1>&2
    ${moduleLoad("bcftools")}
    ulimit -s unlimited
    mkdir -p TMP

    set -x

cat << EOF > TMP/jeter.list
${bcfs.join("\n")}
EOF


    SQRT=`awk 'END{X=NR;z=sqrt(X); print (z==int(z)?z:int(z)+1);}' "TMP/jeter.list"`
    split -a 9 --additional-suffix=.list --lines=\${SQRT} "TMP/jeter.list" TMP/chunck.

    i=1

    find TMP/ -type f -name "chunck.*.list" | while read F
    do
	    bcftools merge -m id -O b -o "TMP/merged.gt.\${i}.bcf" --file-list "\${F}"
            bcftools index --csi "TMP/merged.gt.\${i}.bcf"
	    echo "TMP/merged.gt.\${i}.bcf" >> TMP/to.merge.list
	    i=\$((i+1))
    done



    bcftools merge -m id -O b -o merged.gt.bcf --file-list TMP/to.merge.list
    bcftools index --csi merged.gt.bcf 


	#######################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">merge genotype delly</entry>
		<entry key="versions">${getVersionCmd("bcftools")}</entry>
		<entry key="number of files">${bcfs.size()}</entry>
	</properties>
	EOF

    """

    }

process FILTER_DELLY {
    tag "N=${cases_ctrl_list.size()}"
    cache 'lenient'
    label "process_high"
    memory "5g"
    input:
        val(meta)
        path(delly)
	val(cases_ctrl_list)
	val(merged)
    output:
	path("${params.prefix?:""}sv.bcf"),emit:output
	path("version.xml"),emit:version
	path("${params.prefix?:""}sv.bcf.csi"),emit:index
    script:
	def prefix = params.prefix?:""
	def tag = params.delly2.version.equals("1.1.6")?"":"-t"
    """
    export LC_ALL=C
    ${moduleLoad("bcftools/0.0.0")}
    export PATH=\${PWD}:\${PATH}

    delly filter ${tag} -f germline  -o jeter.bcf "${merged}" 1>&2

    bcftools sort --max-mem "${task.memory.giga}G" -T . -O v -o "jeter1.vcf" jeter.bcf

cat << EOF > jeter.cases.txt
${cases_ctrl_list.findAll{T->T[2].equals("case")}.collect{T->T[1]}.join("\n")}
EOF

cat << EOF > jeter.ctrls.txt
${cases_ctrl_list.findAll{T->T[2].equals("control")}.collect{T->T[1]}.join("\n")}
EOF

    if [  -s "jeter.cases.txt" ] && [  -s "jeter.ctrls.txt"	] ; then
	# rajoute mais pas teste
	bcftools +contrast \
		-0 jeter.ctrls.txt \
		-1 jeter.cases.txt \
		-a PASSOC,FASSOC,NASSOC,NOVELAL,NOVELGT -O v -o jeter2.vcf jeter1.vcf
	mv jeter2.vcf jeter1.vcf
    fi
    rm jeter.cases.txt jeter.ctrls.txt
    

    bcftools  +fill-tags -O v  -o jeter2.vcf jeter1.vcf  -- -t AN,AC,AF,AC_Hom,AC_Het,AC_Hemi,NS
    mv jeter2.vcf jeter1.vcf


    bcftools view -O b -o "${prefix}sv.bcf" jeter1.vcf
    bcftools index "${prefix}sv.bcf"

    rm jeter.bcf jeter1.vcf

	#######################################################################
	cat << EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">delly filter</entry>
		<entry key="versions">${getVersionCmd("bcftools delly")}</entry>
	</properties>
	EOF

    """
    }


process CALL_CNV {
    tag "${name}"
    cache 'lenient'
    afterScript "rm -rf TMP"
    errorStrategy 'finish'
    memory '10 g'
    input:
	val(meta)
	val(genomeId)
	path(delly)
	val(mappability)
	path(bed) //or NO_FILE
	tuple val(name),val(bam),val(sv)
    output:
    	tuple val(name),val(bam),path("${name}.cnv.bcf"),emit:output
	path("version.xml"),emit:version
    script:
        def genome = params.genomes[genomeId]
        def reference =	genome.fasta

	"""
	hostname 1>&2
	mkdir -p TMP
	export TMPDIR=\${PWD}/TMP
	export PATH=\${PWD}:\${PATH}

	delly cnv \
		--outfile "${name}.cnv.bcf" \
		--mappability "${mappability}" \
		--genome "${reference}" \
		--svfile "${sv}" \
		${bed.name.equals("NO_FILE")?"":"--bed-intervals ${bed}"} \\
		"${bam}" 1>&2


	#######################################################################
	cat << EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">call CNV</entry>
		<entry key="sample">${name}</entry>
		<entry key="bam">${bam}</entry>
		<entry key="sv">${sv}</entry>
		<entry key="versions">${getVersionCmd("delly")}</entry>
	</properties>
	EOF
	"""
    }


process MERGE_CNV {
    tag "N=${bcfs.size()}"
    cache 'lenient'
    memory '20 g'
    input:
	val(meta)
	path(delly)
	val(bcfs)
    output:
	path("merged.bcf"),emit:output
	path("version.xml"),emit:version
    script:
    """
    hostname 1>&2
    export LC_ALL=C
    export PATH=\${PWD}:\${PATH}

cat << EOF > jeter.tsv
${bcfs.join("\n")}
EOF

    delly merge -e -p -o merged.bcf  -m 1000 -n 10000000 jeter.tsv 1>&2

rm jeter.tsv

	#######################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">merge CNV</entry>
		<entry key="count">${bcfs.size()}</entry>
		<entry key="version">${getVersionCmd("delly")}</entry>
	</properties>
	EOF

    """
    }


process GENOTYPE_CNV {
    tag "${name}"
    cache 'lenient'
    errorStrategy 'finish'
    afterScript 'rm -f jeter.bcf jeter.bcf.csi'
    memory '10 g'
    input:
	val(meta)
        val(genomeId)
        path(delly)
	val(merged)
        val(mappability)
	path(bed) //or NO_FILE
        tuple val(name),val(bam)
    output:
        path("genotyped.${name}.cnv.bcf"),emit:output
	path("version.xml"),emit:version
    script:
        def genome = params.genomes[genomeId]
        def reference =	genome.fasta
    """
    hostname 1>&2
    mkdir -p TMP
    export TMPDIR=\${PWD}/TMP
    export PATH=\${PWD}:\${PATH}

    delly cnv  \
		--segmentation \
                --vcffile "${merged}" \
		--outfile "jeter.bcf" \
		--mappability "${mappability}" \
		${bed.name.equals("NO_FILE")?"":"--bed-intervals ${bed}"} \\
		--genome "${reference}" \\
		${bam} 1>&2

    mv -v jeter.bcf "genotyped.${name}.cnv.bcf"
    mv -v jeter.bcf.csi "genotyped.${name}.cnv.bcf.csi"
    rm -rf TMP


	#######################################################################
	cat << EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">genotype CNV</entry>
		<entry key="sample">${name}</entry>
		<entry key="bam">${bam}</entry>
		<entry key="version">${getVersionCmd("delly")}</entry>
	</properties>
	EOF
    """
    }

process MERGE_CNV_GENOTYPED {
    tag "N=${bcfs.size()}"
    cache 'lenient'
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
    ${moduleLoad("bcftools")}

cat << EOF > jeter.list
${bcfs.join("\n")}
EOF

    bcftools merge -m id -O b -o merged.gt.bcf --file-list jeter.list
    bcftools index --csi merged.gt.bcf 

	#######################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">merge CNVs</entry>
		<entry key="versions">${getVersionCmd("bcftools")}</entry>
		<entry key="count">${bcfs.size()}</entry>
	</properties>
	EOF


    """
    }

process CLASSIFY_CNV {
    afterScript "rm -f jeter.bcf jeter1.vcf jeter2.vcf"
    cpus 5
    memory "5g"
    input:
	val(meta)
	path(delly)
	val(merged)
    output:
	path("${params.prefix?:""}cnv.bcf"),emit:output
	path("${params.prefix?:""}cnv.bcf.csi"),emit:index
	path("version.xml"),emit:version
    script:
	prefix = params.prefix?:""
    """
    hostname 1>&2
    export LC_ALL=C
    ${moduleLoad("bcftools")}
    export PATH=\${PWD}:\${PATH}

    delly classify -f germline  -o jeter.bcf "${merged}" 1>&2

    bcftools sort --max-mem ${task.memory.giga}G -T . -O b -o "${prefix}cnv.bcf" jeter.bcf
    bcftools index "${prefix}cnv.bcf"

	#######################################################################
	cat << EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">classify CNV</entry>
		<entry key="versions">${getVersionCmd("bcftools delly")}</entry>
	</properties>
	EOF
    """
    }
