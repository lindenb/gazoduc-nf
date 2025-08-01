include {VCF_TO_BED as VCF2BED1    } from '../../modules/bcftools/vcf2bed/main.nf'
include {VCF_TO_BED as VCF2BED2    } from '../../modules/bcftools/vcf2bed/main.nf'

String normContig(String s) {
    if(!s.matches("(chr)?[0-9]+")) return "";
    if(s.startsWith("chr")) s=s.substring(3);
    return s;
    }

workflow PIHAT {
    take:
        meta
        fasta
        fai
        dict
        vcf1kg //[meta,vcf,vcfidx]
        vcfs
    main:


        versions = Channel.empty()
        
        VCF2BED1(vcfs)
        versions = versions.mix(VCF2BED1.out.versions)
        VCF2BED2(vcf1kg)
        versions = versions.mix(VCF2BED2.out.versions)

        ch1 = VCF2BED1.out.output
            .splitCsv(sep:'\t',header:false,elem:1)
            .map{[it[1][0],it[2],it[3]]} /* contig, vcf, vcfidx */
            .map{[normContig(it[0]),it[0],it[1],it[2]]} /* norm_contig, contig, vcf, vcfidx */
            .filter{!it[0].isEmpty()}
            
        ch2 = VCF2BED2.out.output
            .splitCsv(sep:'\t',header:false,elem:1)
            .map{[it[1][0],it[2],it[3]]} /* contig, vcf, vcfidx */
            .map{[normContig(it[0]),it[0],it[1],it[2]]}/* norm_contig, contig, vcf, vcfidx */
            .filter{!it[0].isEmpty()}
            

       

    
        PER_CONTIG(
            fasta,
            fai,
            dict,
            ch1.join(ch2) /* norm_contig, contig, vcf, vcfidx , ikg_contig, 1kgvcf, 1kgidx */
            )
        versions = versions.mix(PER_CONTIG.out.versions)
        
        MERGE(
            fasta,
            fai,
            dict,
            PER_CONTIG.out.bfile
                .collect()
                .map{[[id:"pihat"],it]}
            )
        versions = versions.mix(MERGE.out.versions)


        components = Channel.of(["C1","C2"],["C1","C3"],["C2","C3"])
        formats = Channel.of("pdf","png")
       
        PLOT_MDS(MERGE.out.mds.combine(components).combine(formats))
        versions = versions.mix(PLOT_MDS.out.versions)
       

    emit:
        versions
}


process PER_CONTIG {
tag "${norm_contig}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(norm_contig),
            val(contig_user),path(vcf_user),path(idx_user),
            val(contig_1k),path(vcf_1k),path(idx_1k)
output:
    path("indep_*"),emit:bfile
    path("versions.yml"),emit:versions
script:
    def arg1 = ""
    def half_cpus = task.cpus>2?"--thread ${(task.cpus/2) as int}":""
    def plink_args  = "--const-fid 1 --allow-extra-chr --allow-no-sex --threads ${task.cpus}"
"""
mkdir -p TMP
set -x

##
## ECRACT DATA FROM USER VCF
##
bcftools view  ${half_cpus} --types snps --apply-filters 'PASS,.' --regions "${contig_user}" -m2 -M2 -O u "${vcf_user}" |\\
    bcftools annotate \\
        ${half_cpus}  \\
        -x "INFO,FILTER,QUAL,^FORMAT/GT" \\
        -Ob \\
        -o TMP/jeter1.bcf

bcftools index --threads ${task.cpus} -f   TMP/jeter1.bcf


# view data
bcftools index -s  TMP/jeter1.bcf 1>&2
set +o pipefail
bcftools view -G --no-header TMP/jeter1.bcf |head 1>&2
set -o pipefail



##
## ECRACT DATA FROM 1000 GENOMES
##
echo '${contig_1k}\t${contig_user}' > TMP/chroms.txt 

bcftools view ${half_cpus} --types snps --apply-filters 'PASS,.' --regions "${contig_1k}" -m2 -M2 -O u "${vcf_1k}" |\\
    bcftools annotate \\
        ${half_cpus}  \\
        -x "INFO,FILTER,QUAL,^FORMAT/GT" \\
        --rename-chrs TMP/chroms.txt \\
        -O b \\
        -o TMP/jeter2.bcf


# view data
bcftools index --threads ${task.cpus} -f  TMP/jeter2.bcf
bcftools index -s  TMP/jeter2.bcf 1>&2
set +o pipefail
bcftools view -G --no-header TMP/jeter2.bcf |head 1>&2
set -o pipefail

##
## JOIN
##
bcftools isec --threads ${task.cpus} -O b TMP/jeter1.bcf TMP/jeter2.bcf -p TMP/ISEC -n =2 -w1
find TMP/ISEC 1>&2


# view data
bcftools index -s TMP/ISEC/0000.bcf   1>&2
set +o pipefail
bcftools view -G --no-header TMP/ISEC/0000.bcf  |head 1>&2
set -o pipefail


# rename chromosomes to no chr prefix
paste <(echo '${contig_user}') <(echo '${contig_user}' | sed 's/^chr//') > TMP/chroms.txt 

# select variant after join
bcftools annotate \\
        --threads ${task.cpus} \\
        -i 'F_MISSING < 0.01 && N_ALT=1' \\
        --set-id '%CHROM:%POS:%REF:%ALT' \\
        --rename-chrs TMP/chroms.txt \\
        -O z \\
        -o TMP/jeter1.vcf.gz \\
        TMP/ISEC/0000.bcf      

bcftools index --threads ${task.cpus} -f -t TMP/jeter1.vcf.gz


# check
bcftools index -s  TMP/jeter1.vcf.gz 1>&2
set +o pipefail
bcftools view -G --no-header TMP/jeter1.vcf.gz |head 1>&2
set -o pipefail

rm -rvf TMP/ISEC

# Thank you Floriane S for that wonderful code

plink ${plink_args} \\
    --vcf TMP/jeter1.vcf.gz \\
    --make-bed \\
    --out TMP/data_sansXYMT

## Selection of independents SNP
plink ${plink_args} \\
    --bfile TMP/data_sansXYMT \\
    --indep-pairwise 50 10 0.2 \\
    --out TMP/plink

plink ${plink_args} \\
    --bfile TMP/data_sansXYMT \\
    --extract TMP/plink.prune.in \\
    --make-bed \\
    --out TMP/data_prune

plink ${plink_args} \\
    --bfile TMP/data_prune \\
    --ld-window 50 \\
    --ld-window-kb 5000 \\
    --ld-window-r2 0.2 --r2 \\
    --out TMP/ld

awk '{print \$3}' TMP/ld.ld > TMP/SNP_out.txt

plink ${plink_args} \\
    --bfile TMP/data_prune \\
    --exclude TMP/SNP_out.txt \\
    --make-bed \\
    --out TMP/indep_${norm_contig}

mv TMP/indep_${norm_contig}* ./

cat << EOF > versions.yml
${task.process}:
    plink: todo
	bcftools: "\$(bcftools --version | awk '(NR==1){print \$NF}')"
EOF
"""
}

process MERGE {
tag "${meta.id?:""}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta),path("PLINK/*")
output:
    tuple val(meta),path("*.genome"),emit:genome
    tuple val(meta),path("*.mds"),emit:mds
    path("versions.yml"),emit:versions
script:
    def plink_args  = "--const-fid 1 --allow-extra-chr --allow-no-sex --threads ${task.cpus}"
    def prefix = task.ext.prefix?:meta.id?:"pihat"
"""
mkdir -p TMP
set -x
find PLINK/  -name "*.bim" | sed 's/\\.bim\$//' | sort > TMP/jeter.list
test -s TMP/jeter.list

plink \\
    ${plink_args} \\
    --merge-list TMP/jeter.list \\
    --make-bed \\
    --out TMP/jeter1

## IBS matrix
plink --bfile TMP/jeter1 --genome --out TMP/matIBS_data
wc -l TMP/matIBS_data.genome 1>&2
cp TMP/matIBS_data.genome "${prefix}.genome"



## Remove related individuals
awk '(NR==1 || \$10>0.1) {printf("%s%s\t%s\\n",(NR==1?"#":""),\$3 ,\$4);}' TMP/matIBS_data.genome |\\
    LC_ALL=C sort |\\
    uniq  |\\
    sed 's/^#//' > TMP/indiv_pihat_sup0.1.txt
wc -l TMP/indiv_pihat_sup0.1.txt 1>&2

cp TMP/indiv_pihat_sup0.1.txt ./

plink \\
    ${plink_args} \\
    --bfile TMP/jeter1 \\
    --remove TMP/indiv_pihat_sup0.1.txt \\
    --make-bed \\
    --out TMP/data_sansApp

## MDS
plink \\
    ${plink_args} \\
    --bfile TMP/data_sansApp \\
    --read-genome TMP/matIBS_data.genome \\
    --mds-plot 10 \\
    --cluster \\
    --out TMP/strat_data

mv TMP/strat_data.mds ${prefix}.mds

cat << EOF > versions.yml
${task.process}:
    plink: todo
EOF
"""
}


process PLOT_MDS {
tag "${mds.name} ${C1} ${C2} ${format}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
    tuple val(meta),path(mds),val(C1),val(C2),val(format)
output:
    tuple val(meta),path("*.{pdf,png}"),emit:plot
    path("versions.yml"),emit:versions
script:
    def prefix= task.ext.prefix?:meta.id+".mds"
    def R_main = task.ext.R_main?:"PCA"
    def R_sub = "${C1} / ${C2}"
"""
mkdir -p TMP
cat << '__EOF__' > TMP/jeter.R
strat <- read.table(file="${mds}",header=TRUE)

${format}("TMP/jeter.${format}")
plot(strat\$${C1},strat\$${C2},main ="${R_main}",sub="${R_sub}")
dev.off()
__EOF__

R --no-save < TMP/jeter.R

mv TMP/jeter.${format} "${prefix}.mds.${C1}.${C2}.${format}"

cat << EOF > versions.yml
${task.process}:
    R: todo
EOF
"""
}