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
include {k1_signature      } from '../../modules/utils/k1.nf'
include {SCATTER_TO_BED    } from '../gatk/scatterintervals2bed'

workflow DELLY {
take:
	meta
	fasta
	fai
	dict
	bams_ch
main:
        versions_ch = Channel.empty()
        MAPPABILITY(fai)
        versions_ch = versions_ch.mix(MAPPABILITY.out.versions)

        SCATTER_TO_BED(meta, fasta,fai,dict)
        versions_ch = versions_ch.mix(SCATTER_TO_BED.out.versions)
        
        GET_EXCLUDE(fasta,fai,dict)
        versions_ch = versions_ch.mix(GET_EXCLUDE.out.versions)


        MERGE_EXCLUDE(
            fasta,
            fai,
            SCATTER_TO_BED.out.bed.map{it[1]}.mix(GET_EXCLUDE.out.bed).collect()
            )
        versions_ch = versions_ch.mix(MERGE_EXCLUDE.out.versions)


        /* if there is no meta.status, treat everyone as case */
        bams_ch.map{[it[0].containsKey("status") ? it[0] : it[0].plus("status":"case"), it[1], it[2]] }.
                branch{
                    controls : it[0].status && it[0].status.equals("control")
                    cases : true
                }.set{bams_status_ch}


        CALL_DELLY(fasta, fai,  MERGE_EXCLUDE.out.bed , bams_status_ch.cases )
        versions_ch = versions_ch.mix(CALL_DELLY.out.versions)
        call_vcfs = CALL_DELLY.out.output.map{T->[T[0],T[3],T[4]]}

        MERGE_DELLY(
            fasta,
            fai,
            CALL_DELLY.out.output.
                map{T->[T[3],T[4]]}.
                flatten().
                collect().
                map{[[id:"delly"],it]}
            )
        versions_ch = versions_ch.mix(MERGE_DELLY.out.versions)
      

        GENOTYPE_DELLY(
            fasta,
            fai,
            MERGE_DELLY.out.vcf,
            MERGE_EXCLUDE.out.bed ,
            bams_status_ch.cases.mix(bams_status_ch.controls)
            )
        versions_ch = versions_ch.mix(GENOTYPE_DELLY.out.versions)


        MERGE_GENOTYPES(
            GENOTYPE_DELLY.out.output.
                map{T->[T[1],T[2]]}.
                flatten().
                collect().
                map{[[id:"delly"],it]}
            )
        versions_ch = versions_ch.mix(MERGE_GENOTYPES.out.versions)

        FILTER_DELLY(
            fasta,
            fai,
            MERGE_GENOTYPES.out.vcf
            )
        versions_ch = versions_ch.mix(FILTER_DELLY.out.versions)

emit:
    versions = versions_ch
    vcf = FILTER_DELLY.out.vcf
    call_vcfs
}


process MAPPABILITY {
tag "${fai.name}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	tuple val(meta),path(fai)
output:
	path("blacklist.*"),optional:true,emit:output
    path("versions.yml"),emit:versions
script:
    def k1 = k1_signature()
"""
mkdir -p TMP
	
cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg38}\thttps://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz
1:${k1.hg19}\thttps://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh37.dna.primary_assembly.fa.r101.s501.blacklist.gz
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' |\\
    sed 's/^chr//' |\\
    sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv

join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv |\\
    sort |\\
    uniq > TMP/jeter.url

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

cat << END_VERSIONS > versions.yml
"${task.process}":
    url: \${URL}
END_VERSIONS
"""
}



process GET_EXCLUDE {
tag "${fasta.name}"
label "process_single"
afterScript "rm -f jeter.bed jeter2.bed jeter.interval_list"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
output:
	path("*.bed"),optional:true,emit:bed /* exclude2.bed otherwise collision with bed from scatter_bed */
    path("versions.yml"),emit:versions
script:
 def k1 = k1_signature()
"""
hostname 1>&2
set -o pipefail

mkdir -p TMP
	
cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg19}\thttps://raw.githubusercontent.com/hall-lab/speedseq/master/annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed
1:${k1.hg38}\thttps://gist.githubusercontent.com/chapmanb/4c40f961b3ac0a4a22fd/raw/2025f3912a477edc597e61d911bd1044dc943440/sv_repeat_telomere_centromere.bed
EOF


awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sed 's/^chr//' | sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv

join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv | sort | uniq > TMP/jeter.url


url1=`cat  TMP/jeter.url`


if [ ! -z "\${url1}" ] ; then
	wget -O - "\${url1}" |\\
		cut -f 1,2,3|\\
		jvarkit bedrenamechr -R "${fasta}" --column 1 --convert SKIP > exclude2.bed 
fi

cat << END_VERSIONS > versions.yml
"${task.process}":
    url: \${url1}
END_VERSIONS
"""
}

process MERGE_EXCLUDE {
tag "${fasta.name}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	path("BEDS/*")
output:
	tuple val(meta1),path("exclude_merged.bed"),emit:bed
    path("versions.yml"),emit:versions
script:
    def contig_regex = task.ext.contig_regex?:"^(chr)?[0-9XY]+\$"
"""
hostname 1>&2
mkdir -p TMP
awk -F '\t' '(!(\$1 ~ /${contig_regex}/)) {printf("%s\t0\t%s\\n",\$1,\$2);}'  "${fai}" > TMP/jeter2.bed

cut -f1-3 BEDS/*.bed TMP/jeter2.bed | \
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\
	bedtools merge > exclude_merged.bed

test -s exclude_merged.bed

cat << END_VERSIONS > versions.yml
"${task.process}":
    bedtools: \$(bedtools --version | awk '(NR==1) {print \$NF}')
END_VERSIONS
"""
}


process CALL_DELLY {
    tag "${meta.id}"
    label "process_single" // INCREASE MEMORY IN CONFIG !
    afterScript "rm -rf TMP"
    conda "${moduleDir}/../../conda/delly.yml"
    when:
        task.ext.when == null || task.ext.when
    input:
        tuple val(meta1),path(fasta)
        tuple val(meta2),path(fai)
        tuple val(meta3),path(exclude)
        tuple val(meta),path(bam),path(bai)
    output:
    	tuple val(meta),path(bam),path(bai),path("${meta.id}.bcf"),path("${meta.id}.bcf.csi"),emit:output
        path("versions.yml"),emit:versions
    script:
	    def args1 = exclude?"--exclude ${exclude}":""
        def mapq = task.ext.mapq?:1
        def name = meta.id
        def status = meta.status?:"case"
	"""
	mkdir -p TMP
	export TMPDIR=\${PWD}/TMP

	delly call ${args1} \\
        --map-qual ${mapq} \\
		--outfile "TMP/${name}.bcf" \\
		--genome "${fasta}" \\
		"${bam}" 1>&2

	mv -v TMP/${name}.* ./

cat << END_VERSIONS > versions.yml
"${task.process}":
    delly: \$( delly --version 2>&1 |  awk '(NR==1) {print \$NF;}' )
END_VERSIONS
	"""
    }

//  merge many files: https://github.com/dellytools/delly/issues/158
process MERGE_DELLY {
    tag "${fasta.name}"
    label "process_single"
    afterScript "rm -rf TMP"
    conda "${moduleDir}/../../conda/delly.yml"
    input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta),path("VCFS/*")
    output:
	    tuple val(meta),path("*.bcf"),path("*.bcf.csi"),emit:vcf
        path("versions.yml"),emit:versions
    script:
	    def with_bnd = task.ext.with_bnd?:true
        def prefix = task.ext.prefix?:"merged"
    """
    hostname 1>&2
    export LC_ALL=C
    mkdir -p TMP

    # see https://github.com/dellytools/delly/issues/158
    find VCFS/ -name "*.bcf" > TMP/jeter.tsv
    test -s TMP/jeter.tsv    

    delly merge -o TMP/${prefix}.bcf TMP/jeter.tsv 1>&2

	if ${!with_bnd} ; then
        bcftools view -e 'INFO/SVTYPE="BND"' -O b -o  ${prefix}.bcf TMP/${prefix}.bcf
		bcftools index -f ${prefix}.bcf
	else
		mv TMP/${prefix}.bcf ./
		mv TMP/${prefix}.bcf.csi ./
	fi

cat << END_VERSIONS > versions.yml
"${task.process}":
    delly: \$( delly --version 2>&1 |  awk '(NR==1) {print \$NF;}' )
    bcftools: \$(bcftools version | awk '(NR==1) {print \$NF;}')
END_VERSIONS
    """
    }

process GENOTYPE_DELLY {
    label "process_single" // INCREASE MEMORY IN CONFIG !
    tag "${meta.id}"
    afterScript 'rm -rf TMP'
    conda "${moduleDir}/../../conda/delly.yml"

    input:
	    tuple val(meta1),path(fasta)
	    tuple val(meta2),path(fai)
        tuple val(meta3),path(genotype_vcf),path(genotype_vcf_idx)
        tuple val(meta4),path(exclude)
        tuple val(meta ),path(bam),path(bai)
    output:
        tuple val(meta),path("*bcf",arity:"1"),path("*.csi",arity:"1"),emit:output
        path("versions.yml"),emit:versions
    script:
        def name= meta.id;
        def mapq = task.ext.mapq?:1
    """
    hostname 1>&2
    mkdir -p TMP
    export TMPDIR=\${PWD}/TMP

    delly call \\
        --map-qual ${mapq} \\
        --vcffile "${genotype_vcf}" \\
		--exclude "${exclude}" \\
		--outfile "TMP/jeter.bcf" \\
		--genome "${fasta}" \\
		"${bam}" 1>&2

    mv -v TMP/jeter.bcf "${name}.gt.bcf"
    mv -v TMP/jeter.bcf.csi "${name}.gt.bcf.csi"


cat << END_VERSIONS > versions.yml
"${task.process}":
    delly: \$(delly --version 2>&1 |  awk '(NR==1) {print \$NF;}')
END_VERSIONS
"""
}


process MERGE_GENOTYPES {
    tag "${meta.id}"
    label "process_short"
    afterScript "rm -rf TMP"
    conda "${moduleDir}/../../conda/delly.yml"
    input:
	    tuple val(meta),path("VCFS/*")
    output:
	    tuple val(meta),path("*.bcf",arity:"1"),path("*.csi",arity:"1"),emit:vcf
        path("versions.yml"),emit:versions
    script:
        def prefix = task.ext.prefix?:(meta.id?:"merged")
    // note to self , cd gazoduc-nf if not enough memory
    """
    hostname 1>&2
    ulimit -s unlimited || true
    mkdir -p TMP
    find VCFS/ -name "*.bcf" | sort > TMP/jeter.list
    test -s TMP/jeter.list

    bcftools merge --force-single --threads ${task.cpus} -m id -O b9 -o TMP/${prefix}.bcf --file-list TMP/jeter.list
    bcftools index --threads ${task.cpus} --csi TMP/${prefix}.bcf 
    
    mv TMP/${prefix}.bcf ./
    mv TMP/${prefix}.bcf.csi ./

cat << END_VERSIONS > versions.yml
"${task.process}":
    bcftools: \$(bcftools version | awk '(NR==1) {print \$NF;}')
END_VERSIONS
    """
    }

process FILTER_DELLY {
    tag "${meta.id}"
    conda "${moduleDir}/../../conda/delly.yml"
    label "process_single"
    input:
	    tuple val(meta1),path(fasta)
        tuple val(meta2),path(fai)
        tuple val(meta),path("genotypes.bcf"),path("genotypes.bcf.csi")
    output:
	    tuple val(meta),path("*.bcf",arity:"1"),path("*.csi",arity:"1"),emit:vcf
        path("versions.yml"),emit:versions
    script:
        def args1 = task.ext.args1?:"-t -f germline"
        def prefix = task.ext.prefix?:(meta.id?:"filtered")
    """
    export LC_ALL=C
    mkdir -p TMP

    delly filter ${args1}  -o TMP/jeter.bcf "genotypes.bcf" 1>&2

    bcftools sort --max-mem "${task.memory.giga}G" -T TMP -O v -o "TMP/jeter1.vcf" TMP/jeter.bcf

    bcftools  +fill-tags --threads ${task.cpus} -O v  -o TMP/jeter2.vcf TMP/jeter1.vcf  -- -t AN,AC,AF,AC_Hom,AC_Het,AC_Hemi,NS
    mv TMP/jeter2.vcf TMP/jeter1.vcf

    bcftools view --threads ${task.cpus} -O b9 -o "TMP/${prefix}.bcf" TMP/jeter1.vcf
    bcftools index --threads ${task.cpus} -f "TMP/${prefix}.bcf"

    mv TMP/${prefix}.bcf ./
    mv TMP/${prefix}.bcf.csi ./

cat << END_VERSIONS > versions.yml
"${task.process}":
    delly: \$( delly --version 2>&1 |  awk '(NR==1) {print \$NF;}' )
    bcftools: \$(bcftools version | awk '(NR==1) {print \$NF;}')
END_VERSIONS
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

