/*

Copyright (c) 2025 Pierre Lindenbaum

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
nextflow.enable.dsl=2
//include {ANNOTATE_SV_VCF_01} from "../../subworkflows/annotation/annotation.sv.01.nf"
include { assertKeyExistsAndNotEmpty               } from '../../modules/utils/functions.nf'
include { PREPARE_ONE_REFERENCE                    } from '../../subworkflows/samtools/prepare.one.ref'
include { META_TO_PED                              } from '../../subworkflows/pedigree/meta2ped'
include { MULTIQC                                  } from '../../subworkflows/multiqc'
include { META_TO_BAMS                             } from '../../subworkflows/samtools/meta2bams1'
include { DELLY                                    } from '../../subworkflows/delly2'
include { runOnComplete                            } from '../../modules/utils/functions.nf'
include { READ_SAMPLESHEET                         } from '../../subworkflows/nf/read_samplesheet'
include { JVARKIT_VCFFILTERJDK                     } from '../../modules/jvarkit/vcffilterjdk'
include { JVARKIT_VCF_TO_TABLE as VCF_TO_HTML      } from '../../modules/jvarkit/vcf2table'
include { JVARKIT_VCF_TO_TABLE as VCF_TO_TXT       } from '../../modules/jvarkit/vcf2table'
include { ANNOTSV                                  } from '../../subworkflows/annotsv'
include { BCFTOOLS_QUERY                           } from '../../modules/bcftools/query'
include { PLOT_COVERAGE_01                         } from '../../subworkflows/plotdepth'
include { GTF_INPUT                                } from '../../subworkflows/nf/gtf_input'

workflow {
	  	if(params.fasta==null) {
			throw new IllegalArgumentException("undefined --fasta");
			}
		if(params.samplesheet==null) {
			throw new IllegalArgumentException("undefined --samplesheet");
			}


        versions = Channel.empty()
        multiqc  = Channel.empty()

		 def metadata= [
      		id: "delly2"
      		]
  		PREPARE_ONE_REFERENCE(
			metadata,
			Channel.fromPath(params.fasta).map{f->[[id:f.baseName],f]}
			)
  		versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)
  		
		
		READ_SAMPLESHEET(
			[arg_name:"samplesheet"],
			params.samplesheet
			)
		versions = versions.mix(READ_SAMPLESHEET.out.versions)
  		
  		META_TO_BAMS(
  			metadata ,
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
  			READ_SAMPLESHEET.out.samplesheet
  			)
  		versions = versions.mix(META_TO_BAMS.out.versions)
  		
  		bams_ch = META_TO_BAMS.out.bams
  		
  		META_TO_PED(metadata, bams_ch.map{it[0]})
    	versions = versions.mix(META_TO_PED.out.versions)
    	
		exclude_bed = PREPARE_ONE_REFERENCE.out.complement_bed.first()

		DOWNLOAD_EXCLUDE(
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
			PREPARE_ONE_REFERENCE.out.dict,
			exclude_bed
			)
		versions = versions.mix(DOWNLOAD_EXCLUDE.out.versions)

    	DELLY(
			metadata,
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
			PREPARE_ONE_REFERENCE.out.dict,
			DOWNLOAD_EXCLUDE.out.bed,
			bams_ch
			)
		versions = versions.mix(DELLY.out.versions)
    	
		if(params.jvarkit_vcffilter_script!=null) {
			JVARKIT_VCFFILTERJDK(
				Channel.fromPath(params.jvarkit_vcffilter_script).map{f->[[id:f.baseName],f]},
				META_TO_PED.out.pedigree,
				DELLY.out.vcf.map{meta,vcf,idx->[meta,vcf]}
				)
			versions = versions.mix(JVARKIT_VCFFILTERJDK.out.versions)


			VCF_TO_HTML(
				META_TO_PED.out.pedigree,
				JVARKIT_VCFFILTERJDK.out.vcf
				)
			versions = versions.mix(VCF_TO_HTML.out.versions)
			VCF_TO_TXT(
				META_TO_PED.out.pedigree,
				JVARKIT_VCFFILTERJDK.out.vcf
				)
			versions = versions.mix(VCF_TO_TXT.out.versions)


			BCFTOOLS_QUERY( JVARKIT_VCFFILTERJDK.out.vcf )
			versions = versions.mix(BCFTOOLS_QUERY.out.versions)


			PLOT_COVERAGE_01(
				metadata,
				PREPARE_ONE_REFERENCE.out.fasta,
				PREPARE_ONE_REFERENCE.out.fai,
				PREPARE_ONE_REFERENCE.out.dict,
				BCFTOOLS_QUERY.out.output,
				bams_ch
				)


			ANNOTSV(
				metadata,
				PREPARE_ONE_REFERENCE.out.fasta,
				PREPARE_ONE_REFERENCE.out.fai,
				PREPARE_ONE_REFERENCE.out.dict,
				JVARKIT_VCFFILTERJDK.out.vcf
				)
			versions = versions.mix(ANNOTSV.out.versions)
			multiqc = versions.mix(ANNOTSV.out.multiqc)
			}


    	MULTIQC(
            metadata.plus("id":"delly2"),
            META_TO_PED.out.sample2collection,
            versions,
			[[id:"no_mqc_config"],[]],
            multiqc
            )
		}

runOnComplete(workflow)




workflow xxxx{
		def reference = Channel.fromPath([params.fasta,params.fai,params.dict]).collect()
		delly2_ch = DOWNLOAD_DELLY2()
		gaps_ch = SCATTER_TO_BED(reference)


		map_ch = GET_MAPPABILITY(reference)
		exclude_ch = GET_EXCLUDE(reference)

		xmerge_ch = MERGE_EXCLUDE(reference,gaps_ch.output.mix(exclude_ch.output).collect())


		samplesheet_ch = Channel.fromPath(params.samplesheet).
			splitCsv(header:true,sep:'\t').
			map{it.containsKey("status")?it:it.plus("status":"case")}.
			map{it.containsKey("bai")?it:it.plus("bai":(it.bam.endsWith("bam")?it.bam+".bai":it.bam+".crai"))}.
			map{[it.sample,it.bam,it.bai,it.status]}
			
		cases_ch = samplesheet_ch.filter{it[3].equals("case")}
		controls_ch = samplesheet_ch.filter{!it[3].equals("case")}


		exclude2_ch = exclude_ch.output.ifEmpty(file("NO_FILE"))

		/* run SV discovery */
		if( params.genotype_vcf.equals("NO_FILE")) {

			delly_bcf = CALL_DELLY(reference, delly2_ch.output, exclude2_ch , cases_ch )

			merge_delly = MERGE_DELLY(reference , delly2_ch.output, delly_bcf.output.map{T->T[4]}.collect())
	
			call_vcf = merge_delly.output	
			}
		else
			{
			call_vcf = params.genotype_vcf
			}
	
		genotype_ch = GENOTYPE_DELLY(reference , delly2_ch.output, call_vcf, xmerge_ch.output , samplesheet_ch)	
		merge_gt = MERGE_GENOTYPES(genotype_ch.output.collect())
		filter_delly = FILTER_DELLY( reference , delly2_ch.output, samplesheet_ch.map{it[0]+"\t"+it[3]}.collect() , merge_gt.output ) 

		if(params.with_annotation==true) {
			for_annot_ch = filter_delly.output.map{[
				it[0].name.endsWith(".bcf")?it[0]:it[1], /* vcf */
				it[0].name.endsWith(".bcf")?it[1]:it[0]  /* index */
				]}
			ANNOTATE_SV_VCF_01(reference, for_annot_ch)
			}
}

if (false) {
workflow XX {
	main:







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

	}
}


process GET_MAPPABILITY {
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	path(reference)
output:
	path("blacklist.*"),optional:true,emit:output
script:
	def fai = reference.find{it.name.endsWith(".fai")}.first()
"""
mkdir -p TMP
	
cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:248956422	https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz
1:249250621	https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh37.dna.primary_assembly.fa.r101.s501.blacklist.gz
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sed 's/^chr//' | sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv

join -T TMP -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv | sort | uniq > TMP/jeter.url

URL=`cat  TMP/jeter.url`

if test ! -z "\${URL}"
then
	wget -O TMP/blacklist.gz "\${URL}"
	wget -O TMP/blacklist.gz.fai "\${URL}.fai"
	wget -O TMP/blacklist.gz.gzi "\${URL}.gzi"
	mv TMP/blacklist.* ./
else
	touch blacklist.NO_FILE.gz
	touch blacklist.NO_FILE.gz.fai
	touch blacklist.NO_FILE.gz.gzi
fi
"""
}



process DOWNLOAD_EXCLUDE {
label "process_single"
afterScript "rm -f jeter.bed jeter2.bed jeter.interval_list"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta ),path(complement_bed)
output:
	tuple val(meta ),path("*.bed"),optional:true,emit:bed
	path("versions.yml"),emit:versions
script:
	def url_hg19="https://raw.githubusercontent.com/hall-lab/speedseq/master/annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed"
	def url_hg38=" http://cf.10xgenomics.com/supp/genome/GRCh38/sv_blacklist.bed"
	def url0 =(
		meta1.ucsc_name.equals("hg19")
		?url_hg19
		:meta1.ucsc_name.equals("hg38")
			?url_hg38
			:""
		)
	def url= task.ext.url?:url0
	def prefix = task.ext.prefix?:"exclude"
"""
hostname 1>&2
mkdir -p TMP

cp ${complement_bed} TMP/exclude.bed 

if ${!url.isEmpty()}
then

	curl -L "${url}" |\\
		cut -f 1,2,3|\\
		jvarkit bedrenamechr -R "${fasta}" --column 1 --convert SKIP >> TMP/exclude.bed 

fi

sort -T TMP -t '\t' -k1,1 -k2,2n TMP/exclude.bed |\
	bedtools merge > "${prefix}.bed"

touch versions.yml
"""
stub:
	def prefix = task.ext.prefix?:"exclude"
"""
touch versions.yml ${prefix}.bed
"""
}



process CALL_DELLY {
    tag "${name}"
    label "process_short"
    afterScript "rm -rf TMP"
    input:
	path(reference)
	path(delly)
	path(exclude)
	tuple val(name),path(bam),path(bai),val(status)
    output:
    	tuple val(name),path(bam),path(bai),val(status),path("${name}.bcf"),path("${name}.bcf.csi"),emit:output
    script:
	def fasta = reference.find{it.name.endsWith("a")}
	def args1 = exclude.name.equals("NO_FILE")?"":"--exclude ${exclude}"
	"""
	hostname 1>&2
	mkdir -p TMP
	export TMPDIR=\${PWD}/TMP
	export PATH=\${PWD}:\${PATH}

	delly call ${args1} \
		--outfile "TMP/${name}.bcf" \
		--genome "${fasta}" \
		"${bam}" 1>&2

	mv -v TMP/${name}.* ./
	"""
    }

//  merge many files: https://github.com/dellytools/delly/issues/158
process MERGE_DELLY {
    label "process_short"
    afterScript "rm -rf TMP"
    conda "${moduleDir}/../../conda/bioinfo.01.yml"
    input:
	path(reference)
	path(delly)
	path("VCFS/*")
    output:
	path("merged.*"),emit:output
    script:
        def fasta = reference.find{it.name.endsWith("a")}
	def bnd = params.bnd
    """
    hostname 1>&2
    export LC_ALL=C
    mkdir -p TMP
    export TMPDIR=\${PWD}/TMP
    export PATH=\${PWD}:\${PATH}

    # see https://github.com/dellytools/delly/issues/158

    find VCFS -name "*.bcf" > TMP/jeter.tsv
    test -s TMP/jeter.tsv    

    delly merge -o TMP/merged.bcf TMP/jeter.tsv 1>&2

	if ${!bnd} ; then
                bcftools view -e 'INFO/SVTYPE="BND"' -O b -o  merged.bcf TMP/merged.bcf
		bcftools index -f merged.bcf
	else
		mv TMP/merged.bcf ./
		mv TMP/merged.bcf.csi ./
	fi

    """
    }

process GENOTYPE_DELLY {
    label "process_short"
    tag "${name}"
    cache 'lenient'
    errorStrategy 'finish'
    afterScript 'rm -rf TMP'
    input:
	path(reference)
	path(delly)
        path(merged)
	path(exclude)
        tuple val(name),path(bam),path(bai),val(status)
    output:
        path("genotyped.${name}.bcf*"),emit:output
    script:
	def arg1 = merged.find{it.name.endsWith("f")}
	def fasta = reference.find{it.name.endsWith("a")}
    """
    hostname 1>&2
    mkdir -p TMP
    export TMPDIR=\${PWD}/TMP
    export PATH=\${PWD}:\${PATH}

    delly call --vcffile "${arg1}" \
		--exclude "${exclude}" \
		--outfile "TMP/jeter.bcf" \
		--genome "${fasta}" \
		"${bam}" 1>&2

    mv -v TMP/jeter.bcf "genotyped.${name}.bcf"
    mv -v TMP/jeter.bcf.csi "genotyped.${name}.bcf.csi"
    """
    }


process MERGE_GENOTYPES {
    conda "${moduleDir}/../../conda/bioinfo.01.yml"
    afterScript "rm -rf TMP"
    label "process_high"
    input:
	path(bcfs)
    output:
	path("merged.gt.*"),emit:output
    script:

   if((bcfs instanceof List) &&  bcfs.size() < 1000 )
    """
    hostname 1>&2
    ulimit -s unlimited

cat << EOF > jeter.list
${bcfs.findAll{it.name.endsWith("f")}.join("\n")}
EOF


    bcftools merge -m id -O b -o merged.gt.bcf --file-list jeter.list
    bcftools index --csi merged.gt.bcf 
    rm jeter.list
    """
else
    """
    hostname 1>&2
    ulimit -s unlimited
    mkdir -p TMP

    set -x

    find VCF -name "*f" > TMP/jeter.list

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
    """

    }

process FILTER_DELLY {
    conda "${moduleDir}/../../conda/bioinfo.01.yml"
    label "process_short"
    input:
	path(reference)
        path(delly)
	val(cases_ctrl_list)
	path(merged)
    output:
	path("sv.bcf*"),emit:output
    script:
	                
	def version = params.delly2_version?:""
	def tag = version.equals("1.1.6")?"":"-t"
	def vcfin = merged.find{it.name.endsWith("f")}
    """
    export LC_ALL=C
    export PATH=\${PWD}:\${PATH}
    mkdir -p TMP

    delly filter ${tag} -f germline -o TMP/jeter.bcf "${vcfin}" 1>&2

    bcftools sort --max-mem "${task.memory.giga}G" -T TMP -O v -o "TMP/jeter1.vcf" TMP/jeter.bcf
	               
cat << EOF > TMP/jeter.tsv
${cases_ctrl_list.join("\n")}
EOF

awk -F '\t' '(\$2=="case") {print \$1}' TMP/jeter.tsv > TMP/jeter.cases.txt 
awk -F '\t' '(\$2=="control") {print \$1}' TMP/jeter.tsv > TMP/jeter.ctrls.txt 

    if [ -s "TMP/jeter.cases.txt" ] && [ -s "TMP/jeter.ctrls.txt"	] ; then
	# rajoute mais pas teste
	bcftools +contrast \
		-0 TMP/jeter.ctrls.txt \
		-1 TMP/jeter.cases.txt \
		-a PASSOC,FASSOC,NASSOC,NOVELAL,NOVELGT -O v -o TMP/jeter2.vcf TMP/jeter1.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf
    fi

    bcftools  +fill-tags -O v  -o TMP/jeter2.vcf TMP/jeter1.vcf  -- -t AN,AC,AF,AC_Hom,AC_Het,AC_Hemi,NS
    mv TMP/jeter2.vcf TMP/jeter1.vcf


    bcftools view -O b9 -o "TMP/sv.bcf" TMP/jeter1.vcf
    bcftools index "TMP/sv.bcf"

    mv TMP/sv* ./
    """
    }


process CALL_CNV {
    tag "${name}"
    cache 'lenient'
    afterScript "rm -rf TMP"
     label "process_short"
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

