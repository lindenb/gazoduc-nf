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
include {CARDIOPANEL_DOWNLOAD      } from '../../modules/cardiopanel/download'
include {ALPHAMISSENSE_DOWNLOAD    } from '../../modules/alphamissense/download'
include {AVADA_DOWNLOAD            } from '../../modules/avada/download'
include {CLINVAR_DOWNLOAD          } from '../../modules/clinvar/download'
include {REVEL_DOWNLOAD            } from '../../modules/revel/download'
include {TISSUES_DOWNLOAD          } from '../../modules/jensenlab/tissues/download/main.nf'
include {DISEASES_DOWNLOAD         } from '../../modules/jensenlab/diseases/download/main.nf'
include {SNPEFF_DOWNLOAD           } from '../../modules/snpeff/download/main.nf'
include {BHFUCL_DOWNLOAD           } from '../../modules/bhfucl/download/main.nf'
include {RMSK_DOWNLOAD             } from '../../modules/ucsc/rmsk/download/main.nf'
include {VISTA_DOWNLOAD            } from '../../modules/ucsc/vista/download/main.nf'

workflow ANNOT_SNV {
take:
    meta
    fasta
    fai
    dict
    pedigree
    gtf
    gff3
    vcf // meta, vcf, vcfidx,optional_bed
main:
    versions = Channel.empty()

    //add empty bed if vcf 
    vcf = vcf
        .map{
            if(it.size()==4) return it;
            return [it[0],it[1],it[2],[]]; //add optional bed
        }


    CARDIOPANEL_DOWNLOAD( fasta,fai, dict, gtf  )
    versions = versions.mix(CARDIOPANEL_DOWNLOAD.out.versions)

    ALPHAMISSENSE_DOWNLOAD( fasta,fai, dict, [[id:"nobed"],[]]  )
    versions = versions.mix(ALPHAMISSENSE_DOWNLOAD.out.versions)

    if(fasta[0].ucsc_name && fasta[0].ucsc_name.equals("hg19")) {
        AVADA_DOWNLOAD( fasta,fai, dict )
        versions = versions.mix(AVADA_DOWNLOAD.out.versions)
        avada_vcf = AVADA_DOWNLOAD.out.vcf
    } else 
    {
        avada_vcf = [ [id:"no_avada"],[],[]]
    }

    CLINVAR_DOWNLOAD( fasta,fai, dict, [[id:"nobed"],[]]  )
    versions = versions.mix(CLINVAR_DOWNLOAD.out.versions)

    REVEL_DOWNLOAD( fasta,fai, dict)
    versions = versions.mix(REVEL_DOWNLOAD.out.versions)

    TISSUES_DOWNLOAD(fasta,fai, dict,gtf)
    versions = versions.mix(TISSUES_DOWNLOAD.out.versions)
    
    DISEASES_DOWNLOAD(fasta,fai, dict,gtf)
    versions = versions.mix(DISEASES_DOWNLOAD.out.versions)

    SNPEFF_DOWNLOAD(fai)
    versions = versions.mix(SNPEFF_DOWNLOAD.out.versions)

    BHFUCL_DOWNLOAD(fasta,fai, dict,gtf)
    versions = versions.mix(BHFUCL_DOWNLOAD.out.versions)

    RMSK_DOWNLOAD(fasta,fai, dict)
    versions = versions.mix(RMSK_DOWNLOAD.out.versions)

    VISTA_DOWNLOAD(fasta,fai, dict)
    versions = versions.mix(VISTA_DOWNLOAD.out.versions)


    ANNOTATE(
        fasta,
        fai,
        dict,
        pedigree,
        gtf,
        gff3,
        CARDIOPANEL_DOWNLOAD.out.bed,
        CARDIOPANEL_DOWNLOAD.out.bed_extended,
        ALPHAMISSENSE_DOWNLOAD.out.output,
        avada_vcf,
        CLINVAR_DOWNLOAD.out.vcf,
        REVEL_DOWNLOAD.out.output,
        SNPEFF_DOWNLOAD.out.database,
        TISSUES_DOWNLOAD.out.bed,
        DISEASES_DOWNLOAD.out.bed,
        BHFUCL_DOWNLOAD.out.bed,
        BHFUCL_DOWNLOAD.out.bed_extended,
        RMSK_DOWNLOAD.out.bed,
        VISTA_DOWNLOAD.out.bed,
        vcf
    )

emit:
    versions
    vcf = ANNOTATE.out.vcf
}

process ANNOTATE {
errorStrategy "terminate"
tag "${vcf.name} ${optional_bed?optional_bed.name:""}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
    tuple val(meta1 ),path(fasta)
    tuple val(meta2 ),path(fai)
    tuple val(meta3 ),path(dict)
    tuple val(meta4 ),path(optional_pedigree)
    tuple val(meta5 ),path(gtf),path(gtf_tbi)
    tuple val(meta6 ),path(gff3),path(gff3_tbi)
    /* cardiopanel */
    tuple val(meta7 ),path(cardiopanel),path(cardiopanel_tbi),path(cardiopanel_hdr)
    tuple val(meta8 ),path(cardiopanelx),path(cardiopanelx_tbi),path(cardiopanelx_hdr)
    /** alpha missense */
    tuple val(meta9 ),path(alphamissense),path(alphamissense_idx),path(alphamissense_hdr)
     /** alpha avada */
    tuple val(meta10),path(avada),path(avada_idx)
   /** alpha clinvar */
    tuple val(meta10),path(clinvar),path(clinvar_idx)
    /** revel */
    tuple val(meta11 ),path(revel),path(revel_idx),path(revel_hdr)
    /** snpeff */
    tuple val(meta12 ),path(snpeffdir),path(snpeff_dbname)
    /** tissues */
    tuple val(meta13 ),path(tissues),path(tissues_idx),path(tissues_hdr)
    /** diseases */
    tuple val(meta14 ),path(diseases),path(diseases_idx),path(diseases_hdr)
    /** bhfucl */
    tuple val(meta15 ),path(bhfucl),path(bhfucl_idx),path(bhfucl_hdr)
    tuple val(meta16 ),path(bhfuclx),path(bhfuclx_idx),path(bhfuclx_hdr)
    /** rmsk */
    tuple val(meta17 ),path(rmsk),path(rmsk_idx),path(rmsk_hdr)  
    /** rmvistask */
    tuple val(meta18 ),path(vista),path(vista_idx),path(vista_hdr)  


    tuple val(meta  ),path(vcf),path(vcf_idx),path(optional_bed)
output:
    tuple val(meta),path("*.bcf"),path("*.csi"),emit:vcf
    path("versions.yml"),emit:versions
script:
    def view_args1 = task.ext.args1?:""
    def has_bed = optional_bed?true:false
    def bcftools_norm_args = task.ext.bcftools_norm_args?:""
    def set_id = task.ext.set_id?:""
    def snpeff_args = task.ext.snpeff_args?:""
    def prefix = task.ext.prefix?:vcf.baseName+".annot"
    def vcffilerso_accessions = task.ext.vcffilerso_accessions?:""
    def vcffilerso_args = task.ext.vcffilerso_args?:""
    def gnomadvcf = task.ext.gnomadvcf?:""
    def vcf_gnomad_args = task.ext.vcfgnomad_args?:""
    def gnomad_filterjdk = task.ext.gnomad_filterjdk?:""
"""
mkdir -p TMP
set -x

bcftools query -l ${vcf} | sort | uniq > TMP/samples.txt


################################################################################

if  ${optional_pedigree?true:false}
then

    awk '(\$6=="case" || \$6=="affected") {print \$2;}' '${optional_pedigree}' | sort | uniq > TMP/jeter2.txt
    comm -12 TMP/samples.txt  TMP/jeter2.txt > TMP/jeter.cases.txt

    awk '(\$6=="control" || \$6=="unaffected") {print \$2;}' '${optional_pedigree}' | sort | uniq > TMP/jeter2.txt
    comm -12 TMP/samples.txt TMP/jeter2.txt > TMP/jeter.ctrls.txt

fi

################################################################################

bcftools view \\
    --threads ${task.cpus} \\
    ${view_args1}  \\
    ${has_bed?"--regions-file \"${optional_bed}\"":""} \\
    -O b \\
    -o TMP/jeter2.bcf \\
    "${vcf}"

mv  TMP/jeter2.bcf  TMP/jeter1.bcf

bcftools index -f --threads ${task.cpus}  TMP/jeter1.bcf

################################################################################
if ${!bcftools_norm_args.trim().isEmpty()}
then
	bcftools norm \\
        ${bcftools_norm_args} \\
        --fasta-ref '${fasta}'  \\
        -O u  \\
        TMP/jeter1.bcf |\
        bcftools view  \\
            --write-index \\
            -i 'ALT!=\"*\"' \\
            -o TMP/jeter2.bcf \\

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    mv  TMP/jeter2.bcf.csi  TMP/jeter1.bcf.csi
fi

################################################################################

if ${!set_id.trim().isEmpty()}
then
	bcftools annotate \\
        --write-index \\
         --set-id '${set_id}' \\
        -O b \\
        -o TMP/jeter2.bcf \\
        TMP/jeter1.bcf

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    mv  TMP/jeter2.bcf.csi  TMP/jeter1.bcf.csi
fi


################################################################################

if ${alphamissense?true:false}
then

    bcftools annotate \\
        --write-index \\
        --threads ${task.cpus} \\
        -a "${alphamissense}" \\
        -h "${alphamissense_hdr}" \\
        -c "CHROM,POS,REF,ALT,ALPHAMISSENSE_PATHOGENOCITY,ALPHAMISSENSE_CLASS" \\
        --merge-logic "ALPHAMISSENSE_PATHOGENOCITY:max,ALPHAMISSENSE_CLASS:unique" \\
        -O b \\
        -o TMP/jeter2.bcf \\
        TMP/jeter1.bcf

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    mv  TMP/jeter2.bcf.csi  TMP/jeter1.bcf.csi

fi

################################################################################

if ${rmsk}
then

    bcftools annotate \\
            --write-index \\
            --threads ${task.cpus} \\
            --write-index \\
            -a "${rmsk}" \\
            -h "${rmsk_hdr}" \\
            -c "CHROM,FROM,TO,RMSK"  \\
             -O b \\
        -o TMP/jeter2.bcf \\
        TMP/jeter1.bcf

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    mv  TMP/jeter2.bcf.csi  TMP/jeter1.bcf.csi

fi

################################################################################

if ${vista}
then

    bcftools annotate \\
        --write-index \\
        --threads ${task.cpus} \\
        -a "${vista}" \\
        -h "${vista_hdr}" \\
        -c "CHROM,FROM,TO,VISTA" \\
        --merge-logic 'VISTA:unique' 
        -O b \\
        -o TMP/jeter2.bcf \\
        TMP/jeter1.bcf

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    mv  TMP/jeter2.bcf.csi  TMP/jeter1.bcf.csi

fi



################################################################################

if ${gff3?true:false}
then
    bcftools csq \\
        --threads ${task.cpus} \\
        -O b \\
        --force \\
        --local-csq \\
        --ncsq 10000 \\
        --fasta-ref "${fasta}" \\
        --gff-annot "${gff3}" \\
        -o TMP/jeter2.bcf \\
        TMP/jeter1.bcf

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    bcftools index -f --threads ${task.cpus}  TMP/jeter1.bcf
fi


################################################################################

if [ ! -s "TMP/jeter.cases.txt" ] && [ ! -s "TMP/jeter.ctrls.txt"	]
then

	bcftools +contrast \\
		-0 TMP/jeter.ctrls.txt \\
		-1 TMP/jeter.cases.txt \\
		-a PASSOC,FASSOC,NASSOC,NOVELAL,NOVELGT \\
        -O b \\
        -o TMP/jeter2.bcf \\
        TMP/jeter1.bcf

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    bcftools index -f --threads ${task.cpus}  TMP/jeter1.bcf

fi

################################################################################
if ${cardiopanel?true:false}
then

    bcftools annotate \\
        --threads ${task.cpus} \\
        -a "${cardiopanel}" \\
        -h "${cardiopanel_hdr}" \\
        --write-index \\
        -c "CHROM,POS,END,CARDIOPANEL" \\
        --merge-logic 'CARDIOPANEL:unique' \\
        -O b \\
        -o TMP/jeter2.bcf \\
        TMP/jeter1.bcf

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    bcftools index -f --threads ${task.cpus}  TMP/jeter1.bcf

    bcftools annotate \\
        --threads ${task.cpus} \\
        -a "${cardiopanelx}" \\
        -h "${cardiopanelx_hdr}" \\
        --write-index \\
        --keep-sites -e 'INFO/CARDIOPANEL!= ""' \\
        -c "CHROM,POS,END,CARDIOPANEL_NEAR" \\
        --merge-logic 'CARDIOPANEL_NEAR:unique' \\
        -O b \\
        -o TMP/jeter2.bcf \\
        TMP/jeter1.bcf

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    bcftools index -f --threads ${task.cpus}  TMP/jeter1.bcf

fi

################################################################################

if ${bhfucl?true:false}
then

    bcftools annotate \\
        --threads ${task.cpus} \\
        -a "${bhfucl}" \\
        -h "${bhfucl_hdr}" \\
        --write-index \\
        -c "CHROM,POS,END,BHFUCL" \\
        --merge-logic 'BHFUCL:unique' \\
        -O b \\
        -o TMP/jeter2.bcf \\
        TMP/jeter1.bcf
    
    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    mv  TMP/jeter2.bcf.csi  TMP/jeter1.bcf.csi

    bcftools annotate \\
        --threads ${task.cpus} \\
        -a "${bhfuclx}" \\
        -h "${bhfuclx_hdr}" \\
        --write-index \\
        --keep-sites -e 'INFO/BHFUCL!= ""' \\
        -c "CHROM,POS,END,BHFUCL_NEAR" \\
        --merge-logic 'BHFUCL_NEAR:unique' \\
        -O b \\
        -o TMP/jeter2.bcf \\
        TMP/jeter1.bcf

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    mv  TMP/jeter2.bcf.csi  TMP/jeter1.bcf.csi

fi

################################################################################

if ${avada?true:false}
then

    bcftools annotate \\
        --write-index \\
        --threads ${task.cpus} \\
        -a "${avada}" \\
        -c "AVADA_PMID" \\
        --merge-logic 'AVADA_PMID:unique' \\
        -O b \\
        -o TMP/jeter2.bcf \\
        TMP/jeter1.bcf

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    mv  TMP/jeter2.bcf.csi  TMP/jeter1.bcf.csi
fi

################################################################################

if ${clinvar?true:false}
then
bcftools annotate \\
        --write-index \\
        --threads ${task.cpus} \\
        -a "${clinvar}" \\
        -c "CLINVAR_AF_ESP,CLINVAR_AF_EXAC,CLINVAR_AF_TGP,CLINVAR_ALLELEID,CLINVAR_CLNDN,CLINVAR_CLNDNINCL,CLINVAR_CLNDISDB,CLINVAR_CLNDISDBINCL,CLINVAR_CLNHGVS,CLINVAR_CLNREVSTAT,CLINVAR_CLNSIG,CLINVAR_CLNSIGCONF,CLINVAR_CLNSIGINCL,CLINVAR_CLNVC,CLINVAR_CLNVCSO,CLINVAR_CLNVI,CLINVAR_DBVARID,CLINVAR_GENEINFO,CLINVAR_MC,CLINVAR_ONCDN,CLINVAR_ONCDNINCL,CLINVAR_ONCDISDB,CLINVAR_ONCDISDBINCL,CLINVAR_ONC,CLINVAR_ONCINCL,CLINVAR_ONCREVSTAT,CLINVAR_ONCCONF,CLINVAR_ORIGIN,CLINVAR_RS,CLINVAR_SCIDN,CLINVAR_SCIDNINCL,CLINVAR_SCIDISDB,CLINVAR_SCIDISDBINCL,CLINVAR_SCIREVSTAT,CLINVAR_SCI,CLINVAR_SCIINCL" \\
        -O b \\
        -o TMP/jeter2.bcf \\
        TMP/jeter1.bcf

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    mv  TMP/jeter2.bcf.csi  TMP/jeter1.bcf.csi
fi

################################################################################

if ${revel?true:false}
then

    bcftools annotate \\
         --write-index \\
        --threads ${task.cpus} \\
        -a "${revel}" \\
        -h "${revel_hdr}" \\
        -c "CHROM,POS,REF,ALT,REVEL" \\
        -O b \\
        --merge-logic 'REVEL:max' \\
        -o TMP/jeter2.bcf \\
        TMP/jeter1.bcf

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    mv  TMP/jeter2.bcf.csi  TMP/jeter1.bcf.csi
fi
################################################################################

if ${tissues?true:false}
then

    bcftools annotate \\
            --threads ${task.cpus} \\
            -a "${tissues}" \\
            -h "${tissues_hdr}" \\
            --write-index \\
            -c "CHROM,POS,END,TISSUES_BRENDA,TISSUES" \\
            -O b \\
            --merge-logic 'TISSUES:unique,TISSUES_BRENDA:unique' \\
            -o TMP/jeter2.bcf \\
            TMP/jeter1.bcf

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    mv  TMP/jeter2.bcf.csi  TMP/jeter1.bcf.csi
fi

################################################################################

if ${diseases?true:false}
then

    bcftools annotate \\
            --threads ${task.cpus} \\
            -a "${diseases}" \\
            -h "${diseases_hdr}" \\
            --write-index \\
            -c "CHROM,POS,END,DISEASES_DOID,DISEASES" \\
            -O b \\
            --merge-logic 'DISEASES_DOID:unique,DISEASES:unique' \\
            -o TMP/jeter2.bcf \\
            TMP/jeter1.bcf

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    mv  TMP/jeter2.bcf.csi  TMP/jeter1.bcf.csi
fi




################################################################################

if ${(snpeffdir?true:false) && (snpeff_dbname?true:false)}
then

    bcftools view TMP/jeter1.bcf -O v |\\
        snpEff -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP eff \\
        -dataDir "\${PWD}/${snpeffdir.name}" \\
        -nodownload \\
        ${snpeff_args} \\
        `cat ${snpeff_dbname}` > TMP/jeter1.vcf
    
    bcftools sort --max-mem '${task.memory.giga}G' -T TMP/tmp -O b -o TMP/jeter1.bcf TMP/jeter1.vcf
    bcftools index --threads ${task.cpus} TMP/jeter1.bcf
fi
################################################################################

if ${!vcffilerso_accessions.trim().isEmpty()}
then

    bcftools view TMP/jeter1.bcf -O v |\\
    jvarkit -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP \\
        vcffilterso \\
        ${vcffilerso_args} \\
        --acn "${vcffilerso_accessions}"   >  TMP/jeter1.vcf
    
    bcftools view --threads ${task.cpus} -O b -o TMP/jeter1.bcf  TMP/jeter1.vcf
    bcftools index --threads ${task.cpus} TMP/jeter1.bcf
fi

################################################################################

if ${!gnomadvcf.isEmpty() && !vcf_gnomad_args.isEmpty()}
then

        bcftools view   TMP/jeter1.bcf |\\
        jvarkit -XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP vcfgnomad \\
               ${vcf_gnomad_args} \\
               --gnomad "${gnomadvcf}"  >  TMP/jeter2.vcf
    

    if ${!gnomad_filterjdk.isEmpty()}
    then
        mv  TMP/jeter2.vcf  TMP/jeter1.vcf

        jvarkit -XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP vcffilterjdk \\
                ${gnomad_filterjdk} < TMP/jeter1.vcf > TMP/jeter2.vcf 

    fi

    mv  TMP/jeter2.vcf TMP/jeter1.vcf

    bcftools view --threads ${task.cpus} -O b -o TMP/jeter1.bcf TMP/jeter1.vcf
    bcftools index --threads ${task.cpus} TMP/jeter1.bcf
    rm TMP/jeter1.vcf
fi

################################################################################

mv TMP/jeter1.bcf ${prefix}.bcf
mv TMP/jeter1.bcf.csi ${prefix}.bcf.csi

cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
}