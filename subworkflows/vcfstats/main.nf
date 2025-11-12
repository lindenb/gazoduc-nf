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
include {BCFTOOLS_STATS} from '../../modules/bcftools/stats'

workflow VCF_STATS {
take:
    metadata
    fasta
    fai
    dict
    gtf
    gff3
    pedigree
    sample2pop
    bed
    vcfs
main:
    versions= Channel.empty()
    multiqc = Channel.empty()
    
    if(metadata.with_bcftools_stats==null) {
        log.warn("VCF_STATS: missing metadata.with_bcftools_stats");
        }
    if(metadata.gatk_denovo==null) {
        log.warn("VCF_STATS: missing metadata.gatk_denovo");
        }
    if(metadata.ad_ratio==null) {
        log.warn("VCF_STATS: missing metadata.ad_ratio");
        }
   
    vcf2 = vcfs
        .map{meta,vcf,idx->[meta,[vcf,idx]]}
        .groupTuple()
        .map{meta,files->[meta,files.flatten().sort()]}


    if(metadata.with_bcftools_stats==true) {
        BCFTOOLS_STATS(
            fasta,
            fai,
            bed,
            gtf.map{meta,gff,tbi->[meta,gff]},
            [[id:"no_samples"],[]],
            vcf2
            )
        versions= versions.mix(BCFTOOLS_STATS.out.versions)
        multiqc = multiqc.mix(BCFTOOLS_STATS.out.stats.map{it[1]})
        }

    conditions_ch = Channel.of(
        [
        id: "het",
        label: "heterozygous",
        filter1: "",
        filter2: ".filter(G->G.isHet())"
        ],
        [
        id: "homref",
        label: "homref",
        filter1: "",
        filter2: ".filter(G->G.isHomRef())"
        ],
        [
        id: "homvar",
        label: "homvar",
        filter1: "",
        filter2: ".filter(G->G.isHomVar())"
        ]
    )

    conditions_filter_ch = Channel.of(
        [
        id: "pass",
        label: "PASS",
        filter1: ".filter(V->V.isNotFiltered())",
        filter2: ""
        ],
        [
        id: "filtered",
        label: "FILTERed",
        filter1: ".filter(V->V.isFiltered())",
        filter2: ""
        ]
    )

    if(metadata.gatk_denovo==true) {
        GATK_DE_NOVO(vcf2,bed)
        versions= versions.mix(GATK_DE_NOVO.out.versions)
        multiqc = multiqc.mix(GATK_DE_NOVO.out.multiqc)
        }

    if(metadata.ad_ratio==true) {
        AD_RATIO(vcf2,bed)
        versions= versions.mix(AD_RATIO.out.versions)
        multiqc = multiqc.mix(AD_RATIO.out.multiqc.map{it instanceof List?it:[it]}.flatMap())
        }

    GENOTYPE_QUALITY(vcf2.combine(conditions_ch) ,bed)
    versions= versions.mix(GENOTYPE_QUALITY.out.versions)
    multiqc = multiqc.mix(GENOTYPE_QUALITY.out.multiqc.map{it instanceof List?it:[it]}.flatMap())
    
    FILTERS(vcf2,bed)
    versions= versions.mix(FILTERS.out.versions)
    multiqc = multiqc.mix(FILTERS.out.multiqc)
    

    DEPTH(vcf2.combine(conditions_filter_ch) ,bed)
    versions= versions.mix(DEPTH.out.versions)
    multiqc = multiqc.mix(DEPTH.out.multiqc)


    gatk_field_ch = Channel.of(
        [tag:"FS", range:"0,10,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400"],
        [tag:"MQ", range:"0,10,20,30,40,50,60,70,80,90,100"],
        [tag:"MQRankSum", range:"-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12"],
        [tag:"QD", range:"0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30"],
        [tag:"ReadPosRankSum", range:"-5,-4,-3,-2,-1,0,1,2,3,4,5"],
        [tag:"SOR", range:"0,1,2,3,4,5,6,7,8,9,10"]
    )

    TAG_DISTRIBUTION(vcf2.combine(gatk_field_ch) ,bed)
    versions= versions.mix(TAG_DISTRIBUTION.out.versions)
    multiqc = multiqc.mix(TAG_DISTRIBUTION.out.multiqc)

    GT_FILTERS(vcf2,bed)
    versions= versions.mix(GT_FILTERS.out.versions)
    multiqc = multiqc.mix(GT_FILTERS.out.multiqc)

    SINGLETONS(vcf2,bed)
    versions= versions.mix(SINGLETONS.out.versions)
    multiqc = multiqc.mix(SINGLETONS.out.multiqc)

emit:
    versions
    multiqc
}

process GATK_DE_NOVO {
tag "${meta.id?:""}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"

input:
    tuple val(meta ),path("VCFS/*")
    tuple val(meta2),path(optional_bed)
output:
    path("*_mqc.yml"),optional:true, emit:multiqc
    path("versions.yml"),emit:versions
script:
    def section_id = task.ext.section_id?:"gatk_de_novo_id";
    def section_name = task.ext.section_name?:"VCF stats";
    def prefix  = task.ext.prefix?:"gatkdenovo"
"""
mkdir -p TMP
set -x
set +o pipefail

find VCFS/ \\( -name "*.bcf" -o -name "*.vcf.gz" \\) | sort > TMP/jeter.list


bcftools concat -O u --file-list TMP/jeter.list |\\
    bcftools view --header-only > TMP/header.vcf


if grep -E '##INFO=<ID=(loConfDeNovo|hiConfDeNovo),' TMP/header.vcf
then


cat << EOF > TMP/jeter.code
class SampleInfo {
    final String sn;
    long hiCount=0L;
    long loCount=0L;
    SampleInfo(final String sn) {
        this.sn=sn;
        }
    long total() { return hiCount+loCount;}
    }

final Map<String,SampleInfo> sn2info = new HashMap<>(1000);
EOF

cat << EOF >> TMP/jeter.code
stream().
    forEach(V->{
        for(int side=0;side<2;++side) {
            final String att = (side==0?"loConfDeNovo":"hiConfDeNovo");
            if(!V.hasAttribute(att)) continue;
            for(String sn:V.getAttributeAsStringList(att,"")) {
                SampleInfo si = sn2info.get(sn);
                if(si==null) {
                    si = new SampleInfo(sn);
                    sn2info.put(sn,si);
                    }
                if(side==0) si.loCount++;
                else if(side==1) si.hiCount++;
                }
        }
    });

/* sort sample on number of deNovo, convert to List */
final List<SampleInfo> list = sn2info.values().
    stream().
    sorted((A,B)->Long.compare(B.total(),A.total())).
    collect(Collectors.toList());

for(SampleInfo si:list) {
    println(si.sn+"\t"+si.loCount+"\t"+si.hiCount);
    }

EOF

set -o pipefail

bcftools concat \\
    ${optional_bed?"--regions-file \"${optional_bed}\"":""} \\
    --file-list TMP/jeter.list \\
    --drop-genotypes \\
    -Ov |\\
    jvarkit -Djava.io.tmpdir=TMP bioalcidaejdk -F VCF -f TMP/jeter.code > TMP/jeter.tsv


# just print if there is any DENOVO

if test -s TMP/jeter.tsv
then

cat << EOF > TMP/jeter_mqc.yml
id: "${section_id}"
section_name: "${section_name}"
description: "Number of DeNovo Variant"
plot_type: "bargraph"
pconfig:
    id : "${section_id}_plot"
    title: "GATK Possible DeNovo"
    xlab: "Sample"
    ylab: "Number of Variants"
data:
EOF

awk '{printf("#%s:\\n##loConfDeNovo: %s\\n##hiConfDeNovo: %s\\n",\$1,\$2,\$3);}' TMP/jeter.tsv |\\
    sed 's/#/    /g' >> TMP/jeter_mqc.yml

mv TMP/jeter_mqc.yml ${prefix}_mqc.yml

fi

fi

cat << END_VERSIONS > versions.yml
${task.process}:
    bcftools: \$(bcftools version | awk '(NR==1)  {print \$NF}')
    jvarkit : todo
END_VERSIONS
"""
}

process AD_RATIO {
tag "${meta.id?:""}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
    tuple val(meta ),path("VCFS/*")
    tuple val(meta2),path(optional_bed)
output:
    path("*_mqc.yml"),optional:true, emit:multiqc
    path("versions.yml"),emit:versions
script:
    def section_id = task.ext.section_id?:"ad_ratio_id";
    def section_name = task.ext.section_name?:"VCF stats";
    def prefix  = task.ext.prefix?:"adratio"
"""
mkdir -p TMP
set -x
set +o pipefail

find VCFS/ \\( -name "*.bcf" -o -name "*.vcf.gz" \\) | sort > TMP/jeter.list


bcftools concat -O u --file-list TMP/jeter.list |\\
    bcftools view --header-only > TMP/header.vcf

bcftools query -l TMP/header.vcf > TMP/jeter.samples


if grep -F '##FORMAT=<ID=AD,' TMP/header.vcf && test -s TMP/jeter.samples
then

cat "${moduleDir}/adratio.code" > TMP/jeter.code


bcftools concat \\
    ${optional_bed?"--regions-file \"${optional_bed}\"":""} \\
    --file-list TMP/jeter.list \\
    -Ov |\\
    jvarkit -Djava.io.tmpdir=TMP bioalcidaejdk -F VCF -f TMP/jeter.code > TMP/jeter.tsv



cat << EOF > TMP/template.txt
id: "__TYPE__${section_id}"
section_name: "__TYPE__${section_name}"
description: "AD Ratio for __TYPE__"
plot_type: "bargraph"
pconfig:
    id : "__TYPE__${section_id}_plot"
    title: "AD Ratio __TYPE__"
    xlab: "Sample"
    ylab: "reads with minor allele /DP"
data:
EOF

for F in HET HOM_VAR HOM_REF
do
	sed "s/__TYPE__/\${F}/g"  TMP/template.txt > TMP/jeter_mqc.tsv
	awk -F '\t' -vT=\$F '(\$1==T)' TMP/jeter.tsv | cut -f 2- >> TMP/jeter_mqc.tsv
	
	mv TMP/jeter_mqc.tsv AD_RATIO_\${F}_mqc.yml
done


fi


cat << END_VERSIONS > versions.yml
${task.process}:
    bcftools: \$(bcftools version | awk '(NR==1)  {print \$NF}')
    jvarkit : todo
END_VERSIONS
"""
}


process GENOTYPE_QUALITY {
tag "${meta.id?:""} ${condition.id}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
when:
    task.ext.when == null || task.ext.when
input:
    tuple val(meta ),path("VCFS/*"),val(condition)
    tuple val(meta2),path(optional_bed)
output:
    path("*_mqc.yml"),optional:true, emit:multiqc
    path("versions.yml"),emit:versions
script:
    def section_id = task.ext.section_id?:"gq_id"+condition.id;
    def section_name = task.ext.section_name?:"VCF stats";
    def prefix  = task.ext.prefix?:"gq"
"""
mkdir -p TMP
set -x
set +o pipefail

find VCFS/ \\( -name "*.bcf" -o -name "*.vcf.gz" \\) | sort > TMP/jeter.list


bcftools concat -O u --file-list TMP/jeter.list |\\
    bcftools view --header-only > TMP/header.vcf

bcftools query -l TMP/header.vcf > TMP/jeter.samples


if grep -F '##FORMAT=<ID=GQ,' TMP/header.vcf && test -s TMP/jeter.samples
then

cat "${moduleDir}/gq.code" |\\
    m4 -P \\
        -D__FILTER1__='${condition.filter1}' \\
        -D__FILTER2__='${condition.filter2}' > TMP/jeter.code


bcftools concat \\
    ${optional_bed?"--regions-file \"${optional_bed}\"":""} \\
    --file-list TMP/jeter.list \\
    -Ov |\\
    jvarkit -Djava.io.tmpdir=TMP bioalcidaejdk -F VCF -f TMP/jeter.code > TMP/jeter.tsv


if test -s TMP/jeter.tsv
then

cat << EOF > TMP/jeter_mqc.tsv
id: "${section_id}"
section_name: "${section_name} "
description: "GQ"
plot_type: "bargraph"
pconfig:
    id : "${section_id}_plot"
    title: "Genotype Quality ${condition.label}"
    xlab: "Sample"
    ylab: "Count(Genotype Quality)"
data:
EOF

cat TMP/jeter.tsv >> TMP/jeter_mqc.tsv

mv TMP/jeter_mqc.tsv ${prefix}_mqc.yml


fi

fi


cat << END_VERSIONS > versions.yml
${task.process}:
    bcftools: \$(bcftools version | awk '(NR==1)  {print \$NF}')
    jvarkit : todo
END_VERSIONS
"""
}



process FILTERS {
tag "${meta.id?:""}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
when:
    task.ext.when == null || task.ext.when
input:
    tuple val(meta ),path("VCFS/*")
    tuple val(meta2),path(optional_bed)
output:
    path("*_mqc.yml"),optional:true, emit:multiqc
    path("versions.yml"),emit:versions
script:
    def section_id = task.ext.section_id?:"filters";
    def section_name = task.ext.section_name?:"VCF stats";
    def prefix  = task.ext.prefix?:"filters"
"""
mkdir -p TMP
set -x
set +o pipefail

find VCFS/ \\( -name "*.bcf" -o -name "*.vcf.gz" \\) | sort > TMP/jeter.list


bcftools concat -O u --file-list TMP/jeter.list |\\
    bcftools view --header-only > TMP/header.vcf



if grep -m1 -F '##FILTER=' TMP/header.vcf
then

cat "${moduleDir}/filter.code" |\\
    m4 -P  > TMP/jeter.code


bcftools concat \\
    --drop-genotypes \\
    ${optional_bed?"--regions-file \"${optional_bed}\"":""} \\
    --file-list TMP/jeter.list \\
    -Ov |\\
    jvarkit -Djava.io.tmpdir=TMP bioalcidaejdk -F VCF -f TMP/jeter.code > TMP/jeter.tsv


if test -s TMP/jeter.tsv
then

cat << EOF > TMP/jeter_mqc.tsv
id: "${section_id}"
section_name: "${section_name} "
description: "FILTERS"
plot_type: "bargraph"
pconfig:
    id : "${section_id}_plot"
    title: "Filters
    xlab: "FILTER"
    ylab: "Count"
data:
EOF

awk '{printf("    %s: %s\\n",\$1,\$2);}' TMP/jeter.tsv >> TMP/jeter_mqc.tsv

mv TMP/jeter_mqc.tsv ${prefix}_mqc.yml

fi

fi


cat << END_VERSIONS > versions.yml
${task.process}:
    bcftools: \$(bcftools version | awk '(NR==1)  {print \$NF}')
    jvarkit : todo
END_VERSIONS
"""
}




process DEPTH {
tag "${meta.id?:""} ${condition.id}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
when:
    task.ext.when == null || task.ext.when
input:
    tuple val(meta ),path("VCFS/*"),val(condition)
    tuple val(meta2),path(optional_bed)
output:
    path("*_mqc.yml"),optional:true, emit:multiqc
    path("versions.yml"),emit:versions
script:
    def section_id = task.ext.section_id?:"dp_id"+condition.id;
    def section_name = task.ext.section_name?:"VCF stats";
    def prefix  = task.ext.prefix?:"dp"
"""
mkdir -p TMP
set -x
set +o pipefail

find VCFS/ \\( -name "*.bcf" -o -name "*.vcf.gz" \\) | sort > TMP/jeter.list


bcftools concat -O u --file-list TMP/jeter.list |\\
    bcftools view --header-only > TMP/header.vcf

bcftools query -l TMP/header.vcf > TMP/jeter.samples


if grep -F '##FORMAT=<ID=DP,' TMP/header.vcf && test -s TMP/jeter.samples
then

cat "${moduleDir}/dp.code" |\\
    m4 -P \\
        -D__FILTER1__='${condition.filter1}' \\
        -D__FILTER2__='${condition.filter2}' > TMP/jeter.code


bcftools concat \\
    ${optional_bed?"--regions-file \"${optional_bed}\"":""} \\
    --file-list TMP/jeter.list \\
    -Ov |\\
    jvarkit -Djava.io.tmpdir=TMP bioalcidaejdk -F VCF -f TMP/jeter.code |\\
    sed 's/\\[-Inf/[0/' > TMP/jeter.tsv


if test -s TMP/jeter.tsv
then

cat << EOF > TMP/jeter_mqc.tsv
id: "${section_id}"
section_name: "${section_name} "
description: "DP"
plot_type: "bargraph"
pconfig:
    id : "${section_id}_plot"
    title: "DEPTH ${condition.label}"
    xlab: "Sample"
    ylab: "Count(DP)"
data:
EOF

cat TMP/jeter.tsv >> TMP/jeter_mqc.tsv

mv TMP/jeter_mqc.tsv ${prefix}_mqc.yml


fi

fi


cat << END_VERSIONS > versions.yml
${task.process}:
    bcftools: \$(bcftools version | awk '(NR==1)  {print \$NF}')
    jvarkit : todo
END_VERSIONS
"""
}




process TAG_DISTRIBUTION {
tag "${meta.id?:""} ${condition.tag}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
when:
    task.ext.when == null || task.ext.when
input:
    tuple val(meta ),path("VCFS/*"),val(condition)
    tuple val(meta2),path(optional_bed)
output:
    path("*_mqc.yml"),optional:true, emit:multiqc
    path("versions.yml"),emit:versions
script:
    def section_id = task.ext.section_id?:"dist${condition.tag}";
    def section_name = task.ext.section_name?:"VCF stats";
    def prefix  = task.ext.prefix?:"dist${condition.tag}"
"""
mkdir -p TMP
set -x
set +o pipefail

find VCFS/ \\( -name "*.bcf" -o -name "*.vcf.gz" \\) | sort > TMP/jeter.list


bcftools concat --drop-genotypes -O u --file-list TMP/jeter.list |\\
    bcftools view --header-only > TMP/header.vcf



if grep -m1 -F '##INFO=<ID=${condition.tag},' TMP/header.vcf
then

cat "${moduleDir}/filter.code" |\\
    m4 -P \\
        -D__RANGE__="${condition.range}" \\
        -D__TAG__=${condition.tag} \\
        > TMP/jeter.code


bcftools concat \\
    --drop-genotypes \\
    ${optional_bed?"--regions-file \"${optional_bed}\"":""} \\
    --file-list TMP/jeter.list \\
    -Ov |\\
    jvarkit -Djava.io.tmpdir=TMP bioalcidaejdk -F VCF -f TMP/jeter.code |\\
    sed 's/\\[-Inf/[-/' > TMP/jeter.tsv

if test -s TMP/jeter.tsv
then

cat << EOF > TMP/jeter_mqc.tsv
id: "${section_id}"
section_name: "${section_name} "
description: "${condition.tag}"
plot_type: "bargraph"
pconfig:
    id : "${condition.tag}_plot"
    title: "${condition.tag}"
    xlab: "${condition.tag}"
    ylab: "Count"
data:
EOF

cat TMP/jeter.tsv >> TMP/jeter_mqc.tsv

mv TMP/jeter_mqc.tsv ${prefix}_mqc.yml

fi

fi


cat << END_VERSIONS > versions.yml
${task.process}:
    bcftools: \$(bcftools version | awk '(NR==1)  {print \$NF}')
    jvarkit : todo
END_VERSIONS
"""
}


process GT_FILTERS {
tag "${meta.id?:""}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
when:
    task.ext.when == null || task.ext.when
input:
    tuple val(meta ),path("VCFS/*")
    tuple val(meta2),path(optional_bed)
output:
    path("*_mqc.yml"),optional:true, emit:multiqc
    path("versions.yml"),emit:versions
script:
    def prefix  = task.ext.prefix?:"gtfilter"
"""
mkdir -p TMP
set -x
set +o pipefail

find VCFS/ \\( -name "*.bcf" -o -name "*.vcf.gz" \\) | sort > TMP/jeter.list


bcftools concat -O u --file-list TMP/jeter.list |\\
    bcftools view --header-only > TMP/header.vcf

bcftools query -l TMP/header.vcf > TMP/jeter.samples


if grep -F '##FORMAT=<ID=FT,' TMP/header.vcf && test -s TMP/jeter.samples
then

cat "${moduleDir}/filter.gt.code" |\\
    m4 -P   > TMP/jeter.code


bcftools concat \\
    ${optional_bed?"--regions-file \"${optional_bed}\"":""} \\
    --file-list TMP/jeter.list \\
    -Ov |\\
    jvarkit -Djava.io.tmpdir=TMP bioalcidaejdk -F VCF -f TMP/jeter.code |\\
    sed 's/\\[-Inf/[0/' > TMP/jeter.tsv


if test -s TMP/jeter.tsv
then

mv TMP/jeter.tsv ${prefix}_mqc.yml

fi

fi


cat << END_VERSIONS > versions.yml
${task.process}:
    bcftools: \$(bcftools version | awk '(NR==1)  {print \$NF}')
    jvarkit : todo
END_VERSIONS
"""
}



process SINGLETONS {
tag "${meta.id?:""}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
when:
    task.ext.when == null || task.ext.when
input:
    tuple val(meta ),path("VCFS/*")
    tuple val(meta2),path(optional_bed)
output:
    path("*_mqc.yml"),optional:true, emit:multiqc
    path("versions.yml"),emit:versions
script:
    def prefix  = task.ext.prefix?:"gtfilter"
"""
mkdir -p TMP
set -x
set +o pipefail

find VCFS/ \\( -name "*.bcf" -o -name "*.vcf.gz" \\) | sort > TMP/jeter.list


bcftools concat -O u --file-list TMP/jeter.list |\\
    bcftools view --header-only > TMP/header.vcf

bcftools query -l TMP/header.vcf > TMP/jeter.samples


if grep -F '##FORMAT=<ID=GT,' TMP/header.vcf && test -s TMP/jeter.samples
then

cat "${moduleDir}/sinletons.code" |\\
    m4 -P   > TMP/jeter.code


bcftools concat \\
    ${optional_bed?"--regions-file \"${optional_bed}\"":""} \\
    --file-list TMP/jeter.list \\
    -Ov |\\
    jvarkit -Djava.io.tmpdir=TMP bioalcidaejdk -F VCF -f TMP/jeter.code |\\
    sed 's/\\[-Inf/[0/' > TMP/jeter.tsv

if test -s TMP/jeter.tsv
then

mv TMP/jeter.tsv ${prefix}_mqc.yml

fi

fi


cat << END_VERSIONS > versions.yml
${task.process}:
    bcftools: \$(bcftools version | awk '(NR==1)  {print \$NF}')
    jvarkit : todo
END_VERSIONS
"""
}
