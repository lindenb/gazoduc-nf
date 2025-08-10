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
        sample2pop
        exclude_samples
        exclude_bed
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
            

        no_join= ch1.map{
            [it[0]/* pivot */,it[1]/*ctg*/,it[2]/*vcf*/,it[3]/*idx*/,it[1]/* ctg (same)*/,[]/* no 1kg vcf */,[] /* no 1kg vcf idx */]
            }

        join_ch = ch1.join(ch2) /* norm_contig, contig, vcf, vcfidx , ikg_contig, 1kgvcf, 1kgidx */
            .ifEmpty(no_join) 


        DOWNLOAD_1KG_SAMPLE2POP(meta)
        versions = versions.mix(DOWNLOAD_1KG_SAMPLE2POP.out.versions)


    
        PER_CONTIG(
            fasta,
            fai,
            dict,
            exclude_samples,
            exclude_bed,
            join_ch
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


        MERGE_SAMPLE2POP(
            DOWNLOAD_1KG_SAMPLE2POP.out.tsv,
            sample2pop,
            MERGE.out.merged_plink
            )
        versions = versions.mix(MERGE_SAMPLE2POP.out.versions)


        PLINK_ASSOC(fasta, fai, dict, MERGE.out.merged_plink , MERGE.out.mds )
        versions = versions.mix(PLINK_ASSOC.out.versions)


        PLOT_ASSOC(
            fasta,
            fai,
            dict,
            PLINK_ASSOC.out.assoc
                .flatMap{it[1]}
                .filter{it.name.matches(".*C[123].qassoc")}
                .map{[[id:"pihat"],it]}
            )
        versions = versions.mix(PLOT_ASSOC.out.versions)

        PLOT_PIHAT(MERGE.out.genome)
        versions = versions.mix(PLOT_PIHAT.out.versions)

        components = Channel.of(["C1","C2"],["C1","C3"],["C2","C3"])
        formats = Channel.of("pdf","png")
       
        PLOT_MDS(MERGE.out.mds.combine(components).combine(formats))
        versions = versions.mix(PLOT_MDS.out.versions)
       
        AVERAGE_PIHAT(MERGE.out.genome, MERGE_SAMPLE2POP.out.sample2pop)
        versions = versions.mix(AVERAGE_PIHAT.out.versions)
    emit:
        versions
        genome = MERGE.out.genome
        mds = MERGE.out.mds
}


process PER_CONTIG {
tag "chr${norm_contig}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta4),path(optional_exclude_samples)
    tuple val(meta5),path(optional_exclude_bed)
    tuple val(norm_contig),
            val(contig_user),path(vcf_user),path(idx_user),
            val(optional_contig_1k),path(optional_vcf_1k),path(optional_idx_1k)
output:
    path("indep_*"),emit:bfile
    path("versions.yml"),emit:versions
script:
    def gnomad = task.ext.gnomad?:""
    if(gnomad.isEmpty()) throw new IllegalArgumentException("${task.process} undefined gnomad.");
    def arg1 = ""
    def cpus3 = task.cpus>3?"--thread ${(task.cpus/3) as int}":""
    def plink_args  = "--const-fid 1 --allow-extra-chr --allow-no-sex --threads ${task.cpus}"
"""
mkdir -p TMP
set -x

# 512Mb = max TBI
echo '##INFO=<ID=IN_FILE1,Number=0,Type=Flag,Description="In File 1">' > TMP/file1.header
echo '${contig_user}\t0\t512000000\t1' |bgzip >  TMP/select.bed.gz
tabix -f -p bed TMP/select.bed.gz

##
## ECRACT DATA FROM USER VCF
##
bcftools view  ${cpus3} \\
        --types snps  \\
        --exclude-uncalled \\
        --apply-filters 'PASS,.' \\
        --regions "${contig_user}" \\
        -m2 -M2 \\
        -O u \\
        "${vcf_user}" |\\
    bcftools annotate \\
        ${cpus3}  \\
        -x "INFO,FILTER,QUAL,^FORMAT/GT" \\
        -O u |\\
     bcftools annotate \\
        ${cpus3}  \\
        -a  TMP/select.bed.gz \\
        -h TMP/file1.header \\
        -c "CHROM,FROM,TO,IN_FILE1" \\
        -O b \\
        -o TMP/jeter1.bcf

bcftools index --threads ${task.cpus} -f   TMP/jeter1.bcf


# view data
bcftools index -s  TMP/jeter1.bcf 1>&2
set +o pipefail
bcftools view -G --no-header TMP/jeter1.bcf |head 1>&2
bcftools query -l TMP/jeter1.bcf | cat -n | tail 1>&2
set -o pipefail


if ${optional_vcf_1k?true:false}
then

    ## ADD FLAG to set in source of file 2
    echo '##INFO=<ID=IN_FILE2,Number=0,Type=Flag,Description="In File 2">' > TMP/file2.header


    ##
    ## ECRACT DATA FROM 1000 GENOMES
    ##
    echo '${optional_contig_1k}\t${contig_user}' > TMP/chroms.txt 

    bcftools view ${cpus3} --types snps --apply-filters 'PASS,.' --regions "${optional_contig_1k}" -m2 -M2 -O u "${optional_vcf_1k}" |\\
        bcftools annotate \\
            ${cpus3}  \\
            -x "INFO,FILTER,QUAL,^FORMAT/GT" \\
            --rename-chrs TMP/chroms.txt \\
                -O u |\\
        bcftools annotate \\
            ${cpus3}  \\
            -a  TMP/select.bed.gz \\
            -h TMP/file2.header \\
            -c "CHROM,FROM,TO,IN_FILE2" \\
            -O b \\
            -o TMP/jeter2.bcf
    
    bcftools index --threads ${task.cpus} -f  TMP/jeter2.bcf

    # view data
    bcftools index -s  TMP/jeter2.bcf 1>&2
    set +o pipefail
    bcftools view -G --no-header TMP/jeter2.bcf |head 1>&2
    bcftools query -l TMP/jeter2.bcf | cat -n | tail 1>&2
    set -o pipefail

    ##
    ## MERGE, CHECK VARIANT HAVE BOTH IN_FILE1 and IN_FILE2
    ##
    bcftools merge  -m all -O u \\
        --threads ${task.cpus} \\
        -Ou \\
        TMP/jeter1.bcf TMP/jeter2.bcf |\
        bcftools view \\
            -i 'IN_FILE1=1 && IN_FILE2=1' \\
            -O b -o TMP/jeter3.bcf

    mv TMP/jeter3.bcf TMP/jeter1.bcf

else
        # extract bed for this vcf, extends to avoid too many regions
		bcftools query  \\
            --threads ${task.cpus} \\
            -f '%CHROM\t%POS0\t%END\\n' \\
            TMP/jeter2.bcf |\\
			jvarkit -Djava.io.tmpdir=TMP -jar  bedrenamechr -f "${gnomad}" --column 1 --convert SKIP |\
			sort -T TMP -t '\t' -k1,1 -k2,2n |\
			bedtools merge -i - -d 1000 > TMP/gnomad.bed
	
		if [ ! -s gnomad.bed ] ; then
        		echo "${contig_user}\t0\t1" > TMP/gnomad.bed
		fi


    	# get filtered variants in gnomad
		bcftools query -e 'FILTER=="." || FILTER=="PASS"' \\
            --regions-file TMP/gnomad.bed \\
            -f '%CHROM\t%POS0\t%END\\n' \\
            "${gnomad}" |\\
			jvarkit  -Djava.io.tmpdir=TMP -jar  bedrenamechr -f "${fasta}" --column 1 --convert SKIP |\\
			sort -T	TMP -t '\t' -k1,1	-k2,2n > TMP/x.gnomad.bed
		
		if [ ! -s x.gnomad.bed ] ; then
        		echo "${contig_user}\t0\t1" > TMP/x.gnomad.bed
		fi
		
		bcftools view \\
             --threads ${task.cpus} \\
            --targets-file "^TMP/x.gnomad.bed" \\
            --targets-overlap 2 \\
            -O b \\
            -o TMP/jeter1.bcf \\
            TMP/jeter2.bcf
        
		mv TMP/jeter1.bcf TMP/jeter2.bcf
	
fi

if ${optional_exclude_samples?true:false}
then
    bcftools query -l TMP/jeter1.bcf | sort | uniq > TMP/jeter.a
    sort -T TMP "${optional_exclude_samples}" | uniq > TMP/jeter.b
    # remaining samples
    comm -23 TMP/jeter.a TMP/jeter.b > TMP/jeter.c
    test -s TMP/jeter.c

    # something to remove ?
    if ! cmp TMP/jeter.a TMP/jeter.c
    then
         bcftools view \\
            --threads ${task.cpus} \\
           --samples-file TMP/jeter.c \\
           -O b \\
            -o TMP/jeter2.bcf \\
            TMP/jeter1.bcf
        
        mv TMP/jeter2.bcf TMP/jeter1.bcf
    fi
fi


if ${optional_exclude_bed?true:false}
then
    bcftools view \\
         --threads ${task.cpus} \\
        --targets-file ^${optional_exclude_bed} \\
        --targets-overlap 2 \\
        -O b \\
        -o TMP/jeter2.bcf \\
        TMP/jeter1.bcf
    
    mv TMP/jeter2.bcf TMP/jeter1.bcf
fi

bcftools index --threads ${task.cpus} -f   TMP/jeter1.bcf
 

# rename chromosomes to no chr prefix
paste <(echo '${contig_user}') <(echo '${contig_user}' | sed 's/^chr//') > TMP/chroms.txt 

# select variant after join
bcftools annotate \\
        --threads ${task.cpus} \\
        -i 'F_MISSING < 0.01 && N_ALT=1 && AC>0' \\
        --set-id '%CHROM:%POS:%REF:%ALT' \\
        -x 'INFO' \\
        --rename-chrs TMP/chroms.txt \\
        -O z \\
        -o TMP/jeter1.vcf.gz \\
        TMP/jeter1.bcf


bcftools index --threads ${task.cpus} -f -t TMP/jeter1.vcf.gz


# check
bcftools index -s  TMP/jeter1.vcf.gz 1>&2
set +o pipefail
bcftools view -G --no-header TMP/jeter1.vcf.gz |head 1>&2
bcftools query -l TMP/jeter1.vcf.gz | cat -n | tail 1>&2
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
    tuple val(meta),path("merged.*"),emit:merged_plink
    path("versions.yml"),emit:versions
script:
    def num_mds_components = task.ext.mds_components?:10
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
    --out TMP/merged

## IBS matrix
plink --bfile TMP/merged --genome --out TMP/matIBS_data
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
    --bfile TMP/merged \\
    --remove TMP/indiv_pihat_sup0.1.txt \\
    --make-bed \\
    --out TMP/data_sansApp

## MDS
plink \\
    ${plink_args} \\
    --bfile TMP/data_sansApp \\
    --read-genome TMP/matIBS_data.genome \\
    --mds-plot ${num_mds_components} \\
    --cluster \\
    --out TMP/strat_data

mv TMP/strat_data.mds ${prefix}.mds


mv TMP/merged.* ./

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

process DOWNLOAD_1KG_SAMPLE2POP {
tag "${meta.id?:""}"
afterScript "rm -rf TMP"
label "process_single"
input:
    val(meta)
output:
    tuple val(meta),path("*.tsv"),emit:tsv
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"integrated_call_samples_1kg"
    def url = task.ext.url?:"https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20250704.ALL.ped"
"""
curl -L "${url}" |\\
    cut -f 2,7 |\\
    tail -n+2 >  ${prefix}.tsv

cat << EOF > versions.yml
${task.process}:
    curl: todo
EOF
"""
}

process MERGE_SAMPLE2POP {
tag "${meta1.id?:""}"
afterScript "rm -rf TMP"
label "process_single"
input:
    tuple val(meta1),path(opt_tsv1)
    tuple val(meta2),path(opt_tsv2)
    tuple val(meta),path(merged_plink)
output:
    tuple val(meta1),path("*.tsv"),emit:sample2pop
    path("versions.yml"),emit:versions
script:
    def fam = merged_plink.find{it.name.endsWith(".fam")}
    def prefix = task.ext.prefix?:"merged_sample2pop"
    def other_name = task.ext.other?:"OTHER"
"""
touch jeter.tsv

if ${opt_tsv1?true:false}
then
    cat ${opt_tsv1} >> jeter.tsv
fi

if ${opt_tsv2?true:false}
then
    cat ${opt_tsv2} >> jeter.tsv
fi

cut -f1,2 jeter.tsv |\\
    sort -T . -t '\t' -k1,1 --unique > jeter2.tsv

mv jeter2.tsv jeter.tsv

# all uniq names in 1kg and sample2pop
cut -f1 jeter.tsv | sort -T TMP | uniq >  jeter2.tsv
# sample names in plink.fam
awk '{print \$2}' "${fam}" | sort -T TMP | uniq >  jeter3.tsv
# extract sample without name in jeter.tsv
comm -13 jeter2.tsv jeter3.tsv |\\
    awk '{printf("%s\t${other_name}\\n",\$1);}' >> jeter.tsv

rm jeter2.tsv jeter3.tsv
mv jeter.tsv '${prefix}.tsv'

cat << EOF > versions.yml
${task.process}:
    sort: todo
EOF
"""
}


process PLOT_PIHAT {
tag "${genome.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
    tuple val(meta), path(genome)
output:
    tuple val(meta),path("*.{png,pdf}"),emit:pict
    path("versions.yml"),emit:versions
script:
    def format = task.ext.format?:(meta.format?:"png")
    def prefix=task.ext.prefix?:"pihat"
    def max_pihat = task.ext.max_pihat?:0.1

"""
mkdir -p TMP

cat << '__EOF__' > TMP/jeter.R
genome <- read.table(file="${genome}",header=TRUE)

head(genome)

${format}("TMP/jeter.${format}")
plot(genome\$PI_HAT,
    ylim=c(0,1.0),
    xlab="Individuals Pair",
    ylab="PI-HAT",
    main="PI-HAT"
)
abline(h=${max_pihat},col="blue");
dev.off()
__EOF__

R --no-save < TMP/jeter.R

mv "TMP/jeter.${format}" "${prefix}.samples.${format}"

cat << EOF > versions.yml
${task.process}:
    sort: todo
EOF
"""
}


process AVERAGE_PIHAT {
tag "${genome.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
    tuple val(meta ), path(genome)
    tuple val(meta2), path(sample2group)
output:
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"sample2avg.pihat"
    def sub_title = task.ext.sub?:""
    def max_pihat = task.ext.max_pihat?:0.1
"""
mkdir -p TMP

cat << '__EOF__' > TMP/jeter.awk
(NR>1) {
    P[\$1]+=1.0*(\$10);
    P[\$3]+=1.0*(\$10);
    C[\$1]++;
    C[\$3]++;
    }
END {
    for(S in P) {
        printf("%s\t%f\\n",S,P[S]/C[S]);
        }
    }
__EOF__

# create table sample/avg(pihat)/status
awk -f TMP/jeter.awk '${genome}' |\\
	LC_ALL=C sort -T . -t '\t' -k2,2gr > "TMP/jeter.tsv"


cat << '__EOF__' > TMP/jeter.R
T1<-read.table("TMP/jeter.tsv",sep="\\t",header=FALSE,col.names=c("S","X"),colClasses=c("character","numeric"))

# Read the sample-to-group mapping
sample2group <- read.table("${sample2group}", sep="\t", header=FALSE,col.names=c("S","G"), colClasses=c("character", "character"))

# Merge to get group information for each sample
T1 <- merge(T1, sample2group, by="S", all.x=TRUE)

# Assign a color to each group
groups <- unique(T1\$G)
group_colors <- setNames(rainbow(length(groups)), groups)
T1\$color <- group_colors[T1\$G]

png("TMP/jeter.png")
boxplot(T1\$X ,
    ylim=c(0,max(T1\$X)),
    main="AVG(PIHAT)/SAMPLE",
    sub="${sub_title}",
    xlab="Sample",
    ylab="pihat",
    col=T1\$color
    )
abline(h=${max_pihat},col="blue")


# Add legend
legend("topright", legend=groups, fill=group_colors, title="POP")


dev.off()

__EOF__

R --no-save < TMP/jeter.R

mv TMP/jeter.tsv ${prefix}.tsv

cat << EOF > versions.yml
${task.process}:
    R: todo
EOF
"""
}



process PLINK_ASSOC {
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta4),path("PLINK/*")
	tuple val(meta ),path(mds)
output:
    tuple val(meta ),path("*.qassoc"),emit:assoc
    path("versions.yml"),emit:versions
script:
    def treshold = task.ext.treshold?:5E-8
    def plink_args  = "--const-fid 1 --allow-extra-chr --allow-no-sex --threads ${task.cpus}"


	"""
	hostname 2>&1
	mkdir -p TMP

	# create pheno file https://www.cog-genomics.org/plink/1.9/input#pheno
	awk '(NR==1) {split(\$0,header);} {X=0;for(i=1;i<=NF;i++) {if(header[i] !="SOL") {printf("%s%s",(X==0?"":"\t"),\$i);X=1;}} printf("\\n");}' '${mds}' > TMP/pheno.tsv

	plink \\
        ${plink_args} \\
        --bfile `find PLINK/ -name "*.fam" | sed 's/\\.fam\$//'` \\
		--pheno TMP/pheno.tsv \\
		--all-pheno \\
		--assoc \\
		--out PHENO
	

cat << EOF > versions.yml
${task.process}:
    R: todo
EOF
	"""
	}


process PLOT_ASSOC {
tag "${assoc.name}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
	tuple val(meta ),path(assoc)
output:
    tuple val(meta ),path("*.png"),emit:png
    path("versions.yml"),emit:versions
script:

"""
hostname 1>&2
mkdir -p TMP
set -x

JD1=`which jvarkit`
echo "JD1=\${JD1}" 1>&2
# directory of jvarkit
JD2=`dirname "\${JD1}"`
# find the jar itself
JVARKIT_JAR=`find "\${JD2}/../.." -type f -name "jvarkit.jar" -print -quit`
JVARKIT_DIST=`dirname "\${JVARKIT_JAR}"`


cat "${moduleDir}/Minikit.java" |\\
    sed -e 's/__INPUT__/${assoc}/'  -e 's/__MIN_P_VALUE__/1e-8/' -e 's/__FASTA__/${fasta}/'  > TMP/Minikit.java

javac -d TMP -cp \${JVARKIT_DIST}/jvarkit.jar TMP/Minikit.java
java -Djava.awt.headless=true -cp \${JVARKIT_DIST}/jvarkit.jar:TMP Minikit

#
# Floriane suggested to also add the number of variants per chromosome
#
awk '(\$2!="SNP") {print \$2}' '${assoc}' |\\
	cut -d ':' -f 1 | sort -T TMP | uniq -c |\\
	awk '{printf("%s\t%s\\n",\$2,\$1);}' |\\
	sort -T TMP -t '\t' -k1,1 |\\
	join -t '\t' -1 1 -2 1 - <(sort -t '\t' -k1,1 '${fai}' ) |\\
	awk 'BEGIN{printf("<table class=\\"table\\"><thead><caption>Number of variants per contig</caption><tr><th>CHROM</th><th>Length</th><th>VARIANTS</th><th>Variants per base</th></tr></thead><tbody>\\n");} {printf("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>\\n",\$1,\$3,\$2,\$2/\$3);} END {printf("</tbody></table>\\n");}' >> TMP/jeter.html

mv TMP/jeter.png "${assoc.baseName}.png"

cat << EOF > versions.yml
${task.process}:
    R: todo
EOF
"""
}
