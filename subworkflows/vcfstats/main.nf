include {BCFTOOLS_STATS} from '../../modules/bcftools/stats'

workflow VCF_STATS {
take:
    meta
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
    
   
    vcf2 = vcfs
        .map{[it[0],[it[1],it[2]]]}
        .groupTuple()
        .map{[it[0],it[1].flatten()]}

    BCFTOOLS_STATS(
        fasta,
        fai,
        bed,
        gtf.map{[it[0],it[1]]},
        [[id:"no_samples"],[]],
        vcf2
        )
    versions= versions.mix(BCFTOOLS_STATS.out.versions)
    multiqc = multiqc.mix(BCFTOOLS_STATS.out.stats.map{it[1]})

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


    GATK_DE_NOVO(vcf2,bed)
    versions= versions.mix(GATK_DE_NOVO.out.versions)
    multiqc = multiqc.mix(GATK_DE_NOVO.out.multiqc)

    AD_RATIO(vcf2,bed)
    versions= versions.mix(AD_RATIO.out.versions)
    multiqc = multiqc.mix(AD_RATIO.out.multiqc.map{it instanceof List?it:[it]}.flatMap())


    GENOTYPE_QUALITY(vcf2.combine(conditions_ch) ,bed)
    versions= versions.mix(GENOTYPE_QUALITY.out.versions)
    multiqc = multiqc.mix(GENOTYPE_QUALITY.out.multiqc.map{it instanceof List?it:[it]}.flatMap())
    
     FILTERS(vcf2,bed)
    versions= versions.mix(FILTERS.out.versions)
    multiqc = multiqc.mix(FILTERS.out.multiqc)
    
emit:
    versions
    multiqc
}

process GATK_DE_NOVO {
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
when:
    task.ext.when == null || task.ext.when
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