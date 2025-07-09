include {TRIOS  as TRIO_SNV                       } from '../../subworkflows/trios/main.nf'
include {BCFTOOL_CONCAT                           } from '../../modules/bcftools/concat/main.nf'
include {CLINVAR                                  } from '../../subworkflows/annotation/clinvar/main.nf'
include {ALPHAMISSENSE                            } from '../../subworkflows/annotation/alphamissense/main.nf'
include {BHFUCL                                   } from '../../subworkflows/annotation/bhfucl/main.nf'
include {BCFTOOLS_BCSQ                            } from '../../modules/bcftools/bcsq/main.nf'
include {REVEL                                    } from '../../subworkflows/annotation/revel/main.nf'


workflow WORKFLOW_DENOVO_SNV {
take:
    meta
    fasta
    fai
    dict
    gff3
    gtf
    pedigree
    vcf
main:
    versions = Channel.empty()


    /** get trio data */
    TRIO_SNV(
        meta,
        fasta,
        fai,
        dict,
        pedigree,
        vcf
        )
    vcf = TRIO_SNV.out.vcf
    versions =  versions.mix(TRIO_SNV.out.versions)
    
    KEEP_DE_NOVO(pedigree,vcf)
    versions =  versions.mix(KEEP_DE_NOVO.out.versions)

    BCFTOOL_CONCAT(
        KEEP_DE_NOVO.out.vcf
            .map{[it[1],it[2]]}
            .collect()
            .map{[meta,it.flatten()]},
        [[id:"nobed"],[]
        ])
    versions =  versions.mix(KEEP_DE_NOVO.out.versions)
    vcf = BCFTOOL_CONCAT.out.vcf

    BCFTOOLS_BCSQ(fasta,fai,gff3,vcf )
    vcf = BCFTOOLS_BCSQ.out.vcf
    
    CLINVAR(meta,fasta,fai,dict,[[id:"nobed"],[]],vcf)
    vcf = CLINVAR.out.vcf

    ALPHAMISSENSE(meta,fasta,fai,dict,[[id:"nobed"],[]],vcf)
    vcf = ALPHAMISSENSE.out.vcf

    BHFUCL(meta,fasta,fai,dict,gtf,vcf)
    vcf = BHFUCL.out.vcf
 
    REVEL(meta,fasta,fai,dict,vcf)
    vcf = REVEL.out.vcf

    REPORT(pedigree,vcf)
emit:
    versions
}



process KEEP_DE_NOVO {
    tag "${meta.id}"
    label "process_single"
    conda "${moduleDir}/../../conda/bioinfo.01.yml"
    afterScript "rm -rf TMP"
    input:
        tuple val(met1),path(pedigree)
        tuple val(meta),path(vcf),path(vcfidx)
    output:
        tuple val(meta),path("*.bcf"),path("*.csi"),optional:true,emit:vcf
        path("versions.yml"),emit:versions
    script:
        def prefix = task.ext.prefix?:vcf.baseName+".denovoonly"
    """
    mkdir -p TMP
    
    bcftools view ${vcf} |\\
        jvarkit -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP vcffilterjdk \\
            -e 'return variant.hasAttribute(\\"loConfDeNovo\\")|| variant.hasAttribute(\\"hiConfDeNovo\\") || variant.getAttributeAsInt(\\"MERR\\",0)>0;' |\\
            bcftools view --write-index -O b -o TMP/jeter.bcf

    if test \$(bcftools index -s TMP/jeter.bcf |wc -l) -gt 0
    then
        mv TMP/jeter.bcf ${prefix}.bcf
        mv TMP/jeter.bcf.csi ${prefix}.bcf.csi
    fi

    touch versions.yml
    """
}

process REPORT {
    tag "${meta.id}"
    label "process_single"
    conda "${moduleDir}/../../conda/bioinfo.01.yml"
    afterScript "rm -rf TMP"
    input:
        tuple val(met1),path(pedigree)
        tuple val(meta),path(vcf),path(vcfidx)
    output:
        tuple val(meta),path("*.vcf.gz"),val(meta),path("*.tbi"),emit:vcf
        tuple val(meta),path("*.table.txt")
        tuple val(meta),path("*.genes.tsv")
        tuple val(meta),path("*.bed"),emit:bed //used for igv reports
    script:
        def prefix = task.ext.prefix?:"snv.denovo"
        def args1 = task.ext.args1?:""
    """
    mkdir -p TMP

cat << EOF > TMP/jeter.code
stream().forEach{
    final Set<String> set = new HashSet<>();
    for(int side=0;side<2;side++) {
        final String tag=(i==0?"hiConfDeNovo":"loConfDeNovo");
        if(!variant.hasAttribute(tag)) continue; 
        set.addAll(variant.getAttributeAsStringList(tag,""));
        }
    for(String s : set) {
        println(variant.getContig()+"\t"+(variant.getStart()-1)+"\t"+variant.getEnd()+"\t"+s);
        }
    }
EOF

    bcftools view  ${args1} -Oz -o TMP/jeter.vcf.gz "${vcf}"
    bcftools index -f -t TMP/jeter.vcf.gz 

    bcftools view TMP/jeter.vcf.gz |\\
        jvarkit -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP bioalcidaejdk \\
            -f  TMP/jeter.code |\\
        LC_ALL=C sort -T TMP -k1,1 -k2,2n |\\
        uniq > TMP/jeter.bed
    
    bcftools view TMP/jeter.vcf.gz |\\
        jvarkit -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP vcf2table \\
            --hide 'HOM_REF,NO_CALL' --pedigree ${pedigree} > TMP/jeter.table.txt

    

    bcftools view TMP/jeter.vcf.gz |\\
        jvarkit -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP groupbygene > TMP/jeter.genes.tsv


    mv TMP/jeter.vcf.gz ./${prefix}.vcf.gz
    mv TMP/jeter.vcf.gz.tbi ./${prefix}.vcf.gz.tbi
    mv TMP/jeter.table.txt ./${prefix}.table.txt
    mv TMP/jeter.genes.tsv ./${prefix}.genes.tsv
    mv TMP/jeter.bed  ./${prefix}.bed
    """
}