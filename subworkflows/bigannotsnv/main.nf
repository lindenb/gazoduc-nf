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
include {SIMPLE_REPEATS_DOWNLOAD   } from '../../modules/ucsc/simplerepeats/download/main.nf'
include {VEP_INSTALL_PLUGINS       } from '../../modules/vep/install.plugins'
include {DOWNLOAD_UTR_ANNOTATOR    } from '../../modules/vep/utr.annotator.download'
include {REMAP_DOWNLOAD            } from '../../modules/remap/download'
include {GENCC_DOWNLOAD            } from '../../modules/gencc/download'
include {ENSEMBL_REG_DOWNLOAD      } from '../../modules/ensemblreg/download'
include {HMC_DOWNLOAD              } from '../../modules/hmc/download'
include {GREENDB_DOWNLOAD          } from '../../modules/greendb/download'

String countVariants(def f) {
        return "\necho  \${LINENO} && bcftools query -f '.\\n' \""+f+"\" | wc -l 1>&2" +"\n"
        }


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
    
    if(pedigree[1]) {
        valid_trios = VALID_TRIOS(pedigree).pedigree
       
        versions = versions.mix(VALID_TRIOS.out.versions)
    }  
    else {
        valid_trios =[[id:"no_valid_trio"],[]]
    }
 

    CLINVAR_DOWNLOAD( fasta,fai, dict, [[id:"nobed"],[]]  )
    versions = versions.mix(CLINVAR_DOWNLOAD.out.versions)

    REVEL_DOWNLOAD( fasta,fai, dict)
    versions = versions.mix(REVEL_DOWNLOAD.out.versions)

    TISSUES_DOWNLOAD(fasta,fai, dict,gtf)
    versions = versions.mix(TISSUES_DOWNLOAD.out.versions)
    
    
    DISEASES_DOWNLOAD(fasta,fai, dict,gtf)
    versions = versions.mix(DISEASES_DOWNLOAD.out.versions)

    GENCC_DOWNLOAD(fasta,fai, dict,gtf)
    versions = versions.mix(GENCC_DOWNLOAD.out.versions)

    SNPEFF_DOWNLOAD(fai)
    versions = versions.mix(SNPEFF_DOWNLOAD.out.versions)

    BHFUCL_DOWNLOAD(fasta,fai, dict,gtf)
    versions = versions.mix(BHFUCL_DOWNLOAD.out.versions)

    RMSK_DOWNLOAD(fasta,fai, dict)
    versions = versions.mix(RMSK_DOWNLOAD.out.versions)

    VISTA_DOWNLOAD(fasta,fai, dict)
    versions = versions.mix(VISTA_DOWNLOAD.out.versions)
   
    SIMPLE_REPEATS_DOWNLOAD(fasta,fai, dict)
    versions = versions.mix(SIMPLE_REPEATS_DOWNLOAD.out.versions)

    VEP_INSTALL_PLUGINS(meta)
    versions = versions.mix(VEP_INSTALL_PLUGINS.out.versions)

    DOWNLOAD_UTR_ANNOTATOR(fasta)
    versions = versions.mix(DOWNLOAD_UTR_ANNOTATOR.out.versions)

    REMAP_DOWNLOAD(fasta,fai, dict)
    versions = versions.mix(REMAP_DOWNLOAD.out.versions)

    ENSEMBL_REG_DOWNLOAD(fasta,fai, dict)
    versions = versions.mix(ENSEMBL_REG_DOWNLOAD.out.versions)

    HMC_DOWNLOAD(fasta,fai, dict)
    versions = versions.mix(HMC_DOWNLOAD.out.versions)   

    GREENDB_DOWNLOAD(fasta,fai, dict)
    versions = versions.mix(GREENDB_DOWNLOAD.out.versions)   


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
        VEP_INSTALL_PLUGINS.out.directory,
        valid_trios,
        SIMPLE_REPEATS_DOWNLOAD.out.bed,
        DOWNLOAD_UTR_ANNOTATOR.out.output,
        REMAP_DOWNLOAD.out.bed,
        GENCC_DOWNLOAD.out.bed,
        ENSEMBL_REG_DOWNLOAD.out.bed,
        HMC_DOWNLOAD.out.bed,
        GREENDB_DOWNLOAD.out.bed,
        vcf
    )

    versions = versions.mix(ANNOTATE.out.versions)
emit:
    versions
    vcf = ANNOTATE.out.vcf
}

process ANNOTATE {
errorStrategy "terminate"
tag "${vcf.name} ${optional_bed?optional_bed.name:""}"
label "process_single"
afterScript "rm -rf TMP snpEff_genes.txt snpEff_summary.html"
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
    /** vista */
    tuple val(meta18 ),path(vista),path(vista_idx),path(vista_hdr)  
    /** VEP */
    tuple val(meta19 ),path(vep_plugin_dir) 
    /** valid trios */
    tuple val(meta20 ),path(valid_trios)  
    /** simple repeats */
    tuple val(meta21 ),path(srepeats),path(srepeats_idx),path(srepeats_hdr)  
    /** utr_annotator_file */
    tuple val(meta22 ),path(vep_utr_annotator) 
    /** remap */
    tuple val(meta23 ),path(remap),path(remap_idx),path(remap_hdr)  
    /** gencc */
    tuple val(meta24 ),path(gencc),path(gencc_idx),path(gencc_hdr)
     /** ensemblreg */
    tuple val(meta25 ),path(ensemblreg),path(ensemblreg_idx),path(ensemblreg_hdr)  
    /** HMC */
    tuple val(meta26 ),path(hmc),path(hmc_idx),path(hmc_hdr) 
    /** GREENDB */
    tuple val(meta27 ),path(greendb),path(greendb_idx),path(greendb_hdr)  


    tuple val(meta  ),path(vcf),path(vcf_idx),path(optional_bed)
output:
    tuple val(meta),path("*.bcf"),path("*.csi"),emit:vcf
    path("versions.yml"),emit:versions
script:
    def view_args1 = task.ext.view_args1?:""
    def has_bed = optional_bed?true:false
    def bcftools_norm_args = task.ext.bcftools_norm_args?:""
    def set_id = task.ext.set_id?:""
    def snpeff_args = task.ext.snpeff_args?:""
    def prefix = task.ext.prefix?:vcf.baseName.md5().substring(0,9)+".annot"
    def vcffilerso_accessions = task.ext.vcffilerso_accessions?:""
    def vcffilerso_args = task.ext.vcffilerso_args?:""
    def gnomadvcf = task.ext.gnomadvcf?:""
    def vcf_gnomad_args = task.ext.vcfgnomad_args?:""
    def gnomad_filterjdk = task.ext.gnomad_filterjdk?:""
    def vep_args = task.ext.vep_args?:""
    def vep_spliceai_snv = task.ext.vep_spliceai_snv?:""
    def vep_spliceai_indel = task.ext.vep_spliceai_indel?:""
    def with_vep_spliceai = (task.ext.with_vep_spliceai?true:false) && !vep_spliceai_indel.isEmpty() && !vep_spliceai_snv.isEmpty()
    def vep_loeuf = task.ext.vep_loeuf?:""
    def with_vep_loeuf = (task.ext.with_vep_loeuf?true:false) && !vep_loeuf.isEmpty()
    def vep_gnomad_fields  = task.ext.vep_gnomad_fields?:""
    def with_vep_gnomad = (task.ext.gnomadvcf?true:false) && task.ext.with_vep_gnomad?true:false && !vep_gnomad_fields.isEmpty()
    def cadd = task.ext.cadd?:""
    def with_cadd = (task.ext.with_cadd?true:false) && !cadd.isEmpty()
    def vcfpolyx_size = (task.ext.vcfpolyx_size?:10) 
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
    --write-index \\
    ${view_args1}  \\
    ${has_bed?"--regions-file \"${optional_bed}\"":""} \\
    -O b \\
    -o TMP/jeter2.bcf \\
    "${vcf}"

mv  TMP/jeter2.bcf  TMP/jeter1.bcf
mv  TMP/jeter2.bcf.csi  TMP/jeter1.bcf.csi

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
    ${countVariants("TMP/jeter1.bcf")}
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
    ${countVariants("TMP/jeter1.bcf")}
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
    ${countVariants("TMP/jeter1.bcf")}

fi



################################################################################

if ${hmc?true:false}
then

    bcftools annotate \\
            --write-index \\
            --threads ${task.cpus} \\
            --write-index \\
            -a "${hmc}" \\
            -h "${hmc_hdr}" \\
            -c "CHROM,FROM,TO,HMC"  \\
            --merge-logic 'HMC:min' \\
             -O b \\
            -o TMP/jeter2.bcf \\
        TMP/jeter1.bcf

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    mv  TMP/jeter2.bcf.csi  TMP/jeter1.bcf.csi
    ${countVariants("TMP/jeter1.bcf")}
fi


################################################################################

if ${greendb?true:false}
then

    bcftools annotate \\
            --write-index \\
            --threads ${task.cpus} \\
            --write-index \\
            -a "${greendb}" \\
            -h "${greendb_hdr}" \\
            -c "CHROM,FROM,TO,GREENDB"  \\
            --merge-logic 'GREENDB:unique' \\
             -O b \\
            -o TMP/jeter2.bcf \\
        TMP/jeter1.bcf

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    mv  TMP/jeter2.bcf.csi  TMP/jeter1.bcf.csi
    ${countVariants("TMP/jeter1.bcf")}
fi

################################################################################

if ${rmsk?true:false}
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
    ${countVariants("TMP/jeter1.bcf")}
fi


################################################################################

if ${srepeats?true:false}
then

    bcftools annotate \\
            --write-index \\
            --threads ${task.cpus} \\
            --write-index \\
            -a "${srepeats}" \\
            -h "${srepeats_hdr}" \\
            -c "CHROM,FROM,TO,SREPEAT"  \\
            --merge-logic 'SREPEAT:unique' \\
             -O b \\
            -o TMP/jeter2.bcf \\
        TMP/jeter1.bcf

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    mv  TMP/jeter2.bcf.csi  TMP/jeter1.bcf.csi
    ${countVariants("TMP/jeter1.bcf")}
fi

################################################################################

if ${vista?true:false}
then

    bcftools annotate \\
        --write-index \\
        --threads ${task.cpus} \\
        -a "${vista}" \\
        -h "${vista_hdr}" \\
        -c "CHROM,FROM,TO,VISTA" \\
        --merge-logic 'VISTA:unique' \\
        -O b \\
        -o TMP/jeter2.bcf \\
        TMP/jeter1.bcf

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    mv  TMP/jeter2.bcf.csi  TMP/jeter1.bcf.csi
    ${countVariants("TMP/jeter1.bcf")}
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
    ${countVariants("TMP/jeter1.bcf")}
fi


################################################################################

if [  -s "TMP/jeter.cases.txt" ] && [  -s "TMP/jeter.ctrls.txt"	]
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
    ${countVariants("TMP/jeter1.bcf")}
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
    ${countVariants("TMP/jeter1.bcf")}

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
    ${countVariants("TMP/jeter1.bcf")}
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
    ${countVariants("TMP/jeter1.bcf")}

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
    ${countVariants("TMP/jeter1.bcf")}

fi

################################################################################
if ${valid_trios?true:false} && test -s "${valid_trios}"
then

    awk '{SEX=1;if(\$5=="female") SEX=2;printf("%s\t%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,\$4,SEX);}' "${valid_trios}" > TMP/jeter.trio.ped

    bcftools +mendelian2 \\
        -O z \\
        ${meta1.ucsc_name?(meta1.ucsc_name.equals("hg19")?"--rules GRCh37":meta1.ucsc_name.equals("hg38")?"--rules GRCh38":""):""} \\
        -o TMP/jeter1.vcf.gz \\
        -m a -P TMP/jeter.trio.ped \\
        TMP/jeter1.bcf
    
    bcftools index --threads ${task.cpus} -f -t  TMP/jeter1.vcf.gz
    ${countVariants("TMP/jeter1.vcf.gz")}

    awk -f "${moduleDir}/../../modules/gatk/possibledenovo/pedigree4gatk.awk" "${valid_trios}" > TMP/jeter.trio.ped

    # bug in gatk, it expects that variants are diploids https://github.com/broadinstitute/gatk/blob/342c5ca3adc78e50ab2cc948c71a5ae64574b4ce/src/main/java/org/broadinstitute/hellbender/utils/samples/MendelianViolation.java#L180
    # split haploid, diploid
    bcftools view --threads ${task.cpus} -i 'COUNT(GT="R")==0 && COUNT(GT="A")==0' -O z -o  TMP/jeter1.diploids.vcf.gz TMP/jeter1.vcf.gz
    bcftools index --threads ${task.cpus} -f -t TMP/jeter1.diploids.vcf.gz
    bcftools view --threads ${task.cpus} -e 'COUNT(GT="R")==0 && COUNT(GT="A")==0' -O b -o  TMP/jeter1.haploid.bcf TMP/jeter1.vcf.gz
    bcftools index --threads ${task.cpus} -f  TMP/jeter1.haploid.bcf

    gatk --java-options "-XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" VariantAnnotator \\
        -R "${fasta}" \\
        --annotation PossibleDeNovo \\
        --pedigree  TMP/jeter.trio.ped \\
        -V TMP/jeter1.diploids.vcf.gz \\
        -O "TMP/jeter2.vcf.gz"
    
    bcftools concat -a --threads ${task.cpus} --write-index -O b -o TMP/jeter1.bcf TMP/jeter2.vcf.gz TMP/jeter1.haploid.bcf
    ${countVariants("TMP/jeter1.bcf")}
    rm -f TMP/jeter1.vcf.gz  TMP/jeter1.vcf.gz.tbi TMP/jeter2.vcf.gz TMP/jeter2.vcf.gz.tbi
    rm -f TMP/jeter1.diploids.vcf.gz TMP/jeter1.diploids.vcf.gz.tbi
    rm -f TMP/jeter1.haploid.bcf TMP/jeter1.haploid.bcf.csi
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
    ${countVariants("TMP/jeter1.bcf")}
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
    ${countVariants("TMP/jeter1.bcf")}
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
    ${countVariants("TMP/jeter1.bcf")}
fi
################################################################################

if ${remap?true:false}
then
    bcftools annotate \\
        --write-index \\
        --threads ${task.cpus} \\
        -a "${remap}" \\
        -h "${remap_hdr}" \\
        -c "CHROM,FROM,TO,REMAP" \\
        -o TMP/jeter2.bcf \\
        TMP/jeter1.bcf

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    mv  TMP/jeter2.bcf.csi  TMP/jeter1.bcf.csi
    ${countVariants("TMP/jeter1.bcf")}
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
    ${countVariants("TMP/jeter1.bcf")}
fi


################################################################################


if ${gencc?true:false}
then

    bcftools annotate \\
            --threads ${task.cpus} \\
            -a "${gencc}" \\
            -h "${gencc_hdr}" \\
            --write-index \\
            -c "CHROM,POS,END,GENCC_MONDO,GENCC_DISEASE,GENCC_HPO" \\
            -O b \\
            --merge-logic 'GENCC_MONDO:unique,GENCC_DISEASE:unique,GENCC_HPO:unique' \\
            -o TMP/jeter2.bcf \\
            TMP/jeter1.bcf

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    mv  TMP/jeter2.bcf.csi  TMP/jeter1.bcf.csi
    ${countVariants("TMP/jeter1.bcf")}
fi



################################################################################


if ${ensemblreg?true:false}
then

    bcftools annotate \\
            --threads ${task.cpus} \\
            -a "${ensemblreg}" \\
            -h "${ensemblreg_hdr}" \\
            --write-index \\
            -c "CHROM,POS,END,ENSEMBL_REG" \\
            -O b \\
            --merge-logic 'ENSEMBL_REG:unique' \\
            -o TMP/jeter2.bcf \\
            TMP/jeter1.bcf

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    mv  TMP/jeter2.bcf.csi  TMP/jeter1.bcf.csi
    ${countVariants("TMP/jeter1.bcf")}
fi

################################################################################

    bcftools view  TMP/jeter1.bcf |\\
        jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP vcfpolyx \\
                -n '${vcfpolyx_size}' \\
                --reference '${fasta}' \\
                --tag "POLYX"   |\\
        bcftools view --write-index -O b -o TMP/jeter2.bcf
    
    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    mv  TMP/jeter2.bcf.csi  TMP/jeter1.bcf.csi
    ${countVariants("TMP/jeter1.bcf")}

################################################################################

if ${diseases?true:false}
then

    bcftools annotate \\
            --threads ${task.cpus} \\
            -a "${diseases}" \\
            -h "${diseases_hdr}" \\
            --write-index \\
            -c "CHROM,POS,END,DISEASES,DISEASES_DOID" \\
            -O b \\
            --merge-logic 'DISEASES_DOID:unique,DISEASES:unique' \\
            -o TMP/jeter2.bcf \\
            TMP/jeter1.bcf

    mv  TMP/jeter2.bcf  TMP/jeter1.bcf
    mv  TMP/jeter2.bcf.csi  TMP/jeter1.bcf.csi
    ${countVariants("TMP/jeter1.bcf")}
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
    bcftools index  -f --threads ${task.cpus} TMP/jeter1.bcf
    ${countVariants("TMP/jeter1.bcf")}
fi
################################################################################
if ${!vep_args.isEmpty()}
then


    awk 'END{N=NR;if(N==0) N=1;B=5000.0/N;if(B<1000.0) B=1000;printf("%d\\n",int(B));}' TMP/samples.txt > TMP/jeter.buffer_size


    bcftools view TMP/jeter1.bcf |\\
    vep \\
        ${task.cpus>1?"--fork ${task.cpus}":""} \\
        --buffer_size `cat TMP/jeter.buffer_size` \\
        ${vep_args} \\
        --fasta "${fasta.toRealPath()}" \\
        --dir_plugins "${vep_plugin_dir}/" \\
        --offline \\
        --symbol \\
        --format vcf \\
        --force_overwrite \\
        --sift=b \\
        --polyphen=b \\
        -o TMP/jeter1.vcf \\
        --quiet \\
        --vcf \\
        --no_stats \\
        ${with_vep_gnomad?"--custom ${gnomadvcf},gnomADg,vcf,exact,0,${vep_gnomad_fields}":""} \\
        ${with_vep_spliceai?"--plugin SpliceAI,snv=${vep_spliceai_snv},indel=${vep_spliceai_indel}":""} \\
        ${with_vep_loeuf?"--plugin LOEUF,file=${vep_loeuf},match_by=gene":""} \\
        ${vep_utr_annotator?"--plugin UTRAnnotator,file=${vep_utr_annotator}":""}

    
    bcftools view --write-index --threads ${task.cpus} -O b -o TMP/jeter1.bcf  TMP/jeter1.vcf
    rm TMP/jeter1.vcf
    ${countVariants("TMP/jeter1.bcf")}
fi


#############################################################################################################

if ${!vcffilerso_accessions.trim().isEmpty()}
then

    bcftools view TMP/jeter1.bcf -O v |\\
    java -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP -jar \${HOME}/jvarkit.jar \\
        vcffilterso \\
        ${vcffilerso_args} \\
        --acn "${vcffilerso_accessions}"   >  TMP/jeter1.vcf
    
    bcftools view --write-index --threads ${task.cpus} -O b -o TMP/jeter1.bcf  TMP/jeter1.vcf
    ${countVariants("TMP/jeter1.bcf")}
fi

################################################################################

if ${!gnomadvcf.isEmpty() && !vcf_gnomad_args.isEmpty()}
then

        bcftools view   TMP/jeter1.bcf |\\
        jvarkit -XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP vcfgnomad \\
               ${vcf_gnomad_args} \\
               --gnomad "${gnomadvcf}"  >  TMP/jeter2.vcf
        ${countVariants("TMP/jeter2.vcf")}

    if ${!gnomad_filterjdk.isEmpty()}
    then
        mv  TMP/jeter2.vcf  TMP/jeter1.vcf

        jvarkit -XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP vcffilterjdk \\
                ${gnomad_filterjdk} < TMP/jeter1.vcf > TMP/jeter2.vcf 
        ${countVariants("TMP/jeter2.vcf")}
    fi

    mv  TMP/jeter2.vcf TMP/jeter1.vcf

    bcftools view --write-index --threads ${task.cpus} -O b -o TMP/jeter1.bcf TMP/jeter1.vcf
    rm TMP/jeter1.vcf
    ${countVariants("TMP/jeter1.bcf")}
fi

################################################################################

if ${with_cadd}
then

    bcftools view   TMP/jeter1.bcf |\\
        jvarkit -XX:-UsePerfData  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP vcfcadd \\
            --tabix '${cadd}' > TMP/jeter1.vcf
    
    bcftools view --write-index --threads ${task.cpus} -O b -o TMP/jeter1.bcf TMP/jeter1.vcf
    rm TMP/jeter1.vcf
    ${countVariants("TMP/jeter1.bcf")}
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


process VALID_TRIOS {
executor "local"
afterScript "rm -rf TMP"
input:
    tuple val(meta),path(ped)
output:
    tuple val(meta),path("trios.ped.tsv"),emit:pedigree
    path("versions.yml"),emit:versions
script:
"""
mkdir -p TMP
tr " " "\t" < "${ped}" | tr -s "\t" |sort -t '\t' -k2,2 | uniq > TMP/samples.tsv

awk '(\$3!="0" && \$4!="0")'  TMP/samples.tsv  > TMP/trios.tsv
cut -f3 TMP/trios.tsv | sort | uniq > TMP/father.txt
cut -f4 TMP/trios.tsv | sort | uniq > TMP/mother.txt

join -t '\t' -1 2 -2 1 -o '1.1,1.2,1.3,1.4,1.5' TMP/samples.tsv TMP/father.txt >> TMP/trios.tsv
join -t '\t' -1 2 -2 1 -o '1.1,1.2,1.3,1.4,1.5' TMP/samples.tsv TMP/mother.txt >> TMP/trios.tsv

sort -t '\t' -k2,2  TMP/trios.tsv | uniq > trios.ped.tsv

cat << END_VERSIONS > versions.yml
"${task.process}":
	join: todo
    awk: todo
END_VERSIONS
"""
}
