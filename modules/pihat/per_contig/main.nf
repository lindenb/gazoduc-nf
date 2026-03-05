/*

Copyright (c) 2026 Pierre Lindenbaum

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
include { verify       } from '../../../modules/utils/functions.nf'
include { isBlank      } from '../../../modules/utils/functions.nf'
include { parseBoolean } from '../../../modules/utils/functions.nf'

process PER_CONTIG {
tag "${meta.id}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta4),path(optional_exclude_samples)
    tuple val(meta5),path(optional_exclude_bed)
    tuple val(metaK),path(optional_vcf_1k),path(optional_idx_1k) // optional
    tuple val(metaG),path(gnomad),path(gnomad_idx) // optional
    tuple val(meta ),path(vcf_user),path(idx_user)
output:
    tuple val(meta ),path("*.bim"),path("*.bed"),path("*.fam"),emit:bfile
    path("versions.yml"),emit:versions
script:
    def contig_user = meta.contig
     verify(!isBlank(contig_user),"${task.process} contig_user is blank")
    def norm_contig = meta.norm_contig
    verify(!isBlank(norm_contig),"${task.process} norm_contig is blank")
    def optional_contig_1k = metaK.contig?:"${contig_user}"
    def optional_contig_gnomad = metaG.contig?:"${contig_user}"
    def arg1 = ""
    def cpus3 = task.cpus>3?"--thread ${(task.cpus/3) as int}":""
    def plink_args  = "--const-fid 1 --allow-extra-chr --allow-no-sex --threads ${task.cpus}"
    def no_call_to_hom_ref= parseBoolean(task.ext.no_call_to_hom_ref?:true)
    def prefix = task.ext.prefix?:"${meta.id}.indep"
    def jvm = task.ext.jvm?:"-Djava.io.tmpdir=TMP"
    def maf = task.ext.maf?:0.1
    def f_missing = task.ext.f_missing?:0.01
"""
mkdir -p TMP
set -x

function dump_vcf() {
    # view data
    echo "DEBUG VCF for \$1:" 1>&2
    bcftools index -s  "\$1" 1>&2
    bcftools +fill-tags -O u "\$1" -- -t  AN,AC,AF,AC_Hom,AC_Het,AC_Hemi,NS |\\
        bcftools view -G --no-header |head 1>&2 || true
    bcftools query -l "\$1" | cat -n | tail 1>&2 || true
    }

echo "${contig_user}\t${norm_contig}" > TMP/rename_contig.tsv

##
## ECRACT DATA FROM USER VCF
##
bcftools view  ${cpus3} \\
        --types snps  \\
	--min-af "${maf}" --max-af "${1.0 - (maf as double)}"  \
        -i 'F_MISSING < ${f_missing}' \\
        --exclude-uncalled \\
        --apply-filters 'PASS,.' \\
        --regions "${contig_user}" \\
        -m2 -M2 \\
        -O u \\
        "${vcf_user}" |\\
    bcftools view \\
        -O u \\
        ${optional_exclude_bed?"--targets-file ^${optional_exclude_bed} --targets-overlap 2 ":""} |\\
    bcftools annotate \\
        ${cpus3}  \\
        ${contig_user==norm_contig?"":"--rename-chrs TMP/rename_contig.tsv"} \\
        -x "INFO,ID,FILTER,QUAL,^FORMAT/GT" \\
        -O b \\
        -o TMP/jeter1.bcf
       

bcftools index --threads ${task.cpus} -f   TMP/jeter1.bcf

dump_vcf TMP/jeter1.bcf


# 
# 1000 GENOMES DATA
#
if ${optional_vcf_1k?true:false}
then

    bcftools query -f '%CHROM\t%POS0\t%END\\n' TMP/jeter1.bcf > TMP/jeter.bed
    test -s TMP/jeter.bed


    ##
    ## ECRACT DATA FROM 1000 GENOMES
    ##
    echo '${optional_contig_1k}\t${norm_contig}' > TMP/chroms.txt 

    bcftools view ${cpus3} \\
        --types snps \\
        --apply-filters 'PASS,.' \\
        --regions "${optional_contig_1k}" \\
        -m2 -M2 \\
        -O u \\
        "${optional_vcf_1k}" |\\
         bcftools annotate \\
            ${cpus3}  \\
            -x "INFO,ID,FILTER,QUAL,^FORMAT/GT" \\
            --rename-chrs TMP/chroms.txt \\
            -O u |\\
            bcftools view \\
                -O b \\
                --write-index \\
                --targets-file  TMP/jeter.bed --targets-overlap 2 \\
                -o TMP/jeter2.bcf
        

    mkdir -p TMP/ISEC
    bcftools isec --threads ${task.cpus} --write-index -O b  -c none -n=2 -w1,2 -p TMP/ISEC TMP/jeter1.bcf TMP/jeter2.bcf
    find TMP/ISEC 1>&2
    
    ##
    ## MERGE, VARIANT HAVE BOTH IN_FILE1 and IN_FILE2
    ##
    bcftools merge  -m all -O u \\
        --threads ${task.cpus} \\
       -O u \\
        TMP/ISEC/0000.bcf TMP/ISEC/0001.bcf |\\
            bcftools view -i 'F_MISSING<0.01'  -O b -o TMP/jeter3.bcf \\

    rm -rvf TMP/ISEC

    mv TMP/jeter3.bcf TMP/jeter1.bcf
    bcftools index --threads ${task.cpus} -f   TMP/jeter1.bcf


    dump_vcf TMP/jeter1.bcf
fi



#
# REMOVE VARIANTS OF LOW QUALITY IN GNOMAD
#
if ${gnomad?true:false}
then

        # extract bed for this vcf, extends to avoid too many regions
        bcftools query  \\
            -f '%CHROM\t%POS0\t%END\\n' \\
            TMP/jeter1.bcf |\\
			jvarkit ${jvm} bedrenamechr -f "${gnomad}" --column 1 --convert SKIP |\\
			sort -T TMP -t '\t' -k1,1 -k2,2n |\\
			bedtools merge -i - -d 1000 > TMP/gnomad.bed
	
		if [ ! -s gnomad.bed ] ; then
        		echo "${optional_contig_gnomad}\t0\t1" > TMP/gnomad.bed
		fi


    	# get filtered variants in gnomad
		bcftools query -e 'FILTER=="." || FILTER=="PASS"' \\
            --regions-file TMP/gnomad.bed \\
            -f '%CHROM\t%POS0\t%END\\n' \\
            "${gnomad}" |\\
			jvarkit ${jvm} bedrenamechr -f "${fasta}" --column 1 --convert SKIP |\\
			sort -T	TMP -t '\t' -k1,1	-k2,2n > TMP/x.gnomad.bed
		
		if [ ! -s x.gnomad.bed ] ; then
        		echo "${contig_user}\t0\t1" > TMP/x.gnomad.bed
		fi
		
		bcftools view \\
             --threads ${task.cpus} \\
            --targets-file "^TMP/x.gnomad.bed" \\
            --targets-overlap 2 \\
            -O b \\
            -o TMP/jeter2.bcf \\
            TMP/jeter1.bcf
        	
	mv TMP/jeter2.bcf TMP/jeter1.bcf
    bcftools index --threads ${task.cpus} -f   TMP/jeter1.bcf

	dump_vcf TMP/jeter1.bcf
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

#
# translate NO_CALL to HOM_REF
#
if ${no_call_to_hom_ref}
then
    bcftools +setGT --threads ${task.cpus} -O b  -o "TMP/jeter2.bcf"  TMP/jeter1.bcf  --  -t '.' -n 0
    mv TMP/jeter2.bcf TMP/jeter1.bcf
fi

#
# Assign variant ID
#
bcftools annotate --threads ${task.cpus} --set-id +'%VKX' -O z -o TMP/jeter1.vcf.gz TMP/jeter1.bcf
bcftools index --threads ${task.cpus} -f -f   TMP/jeter1.vcf.gz
dump_vcf TMP/jeter1.vcf.gz


# Thank you Floriane S for that wonderful code
# convert VCF to plink (BCF marche pas ?)

plink ${plink_args} \\
	--vcf TMP/jeter1.vcf.gz \\
	--make-bed \\
	--out TMP/jeter1 1>&2

find TMP -type f 1>&2

#
# Dans le fichier hardyweinberg.hwe enlever les variants dont la colonne P < 0.00001
#

plink ${plink_args} \\
	--bfile TMP/jeter1 \\
	--hardy gz \\
	--out "TMP/hardyweinberg.txt" 1>&2

find TMP -type f 1>&2
file TMP/hardyweinberg.txt.hwe.gz 1>&2

gunzip -c TMP/hardyweinberg.txt.hwe.gz |\\
	awk '(\$9 <  0.00001) {print \$2}'  > TMP/xclude_ids.txt

wc -l TMP/xclude_ids.txt 1>&2

# remove variants with those IDS
plink ${plink_args} \\
	--bfile TMP/jeter1 \\
	--make-bed \\
	--exclude TMP/xclude_ids.txt \\
	--out TMP/data_sansXYMT 1>&2


## Selection of independents SNP
plink2 ${plink_args} \\
    --bfile TMP/data_sansXYMT \\
    --indep-pairwise 50 10 0.2 \\
    --out TMP/plink

plink2 ${plink_args} \\
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

plink2 ${plink_args} \\
    --bfile TMP/data_prune \\
    --exclude TMP/SNP_out.txt \\
    --make-bed \\
    --out TMP/${prefix}

mv TMP/${prefix}* ./

cat << EOF > versions.yml
${task.process}:
    plink: "\$(plink --version | awk '{print \$2}')"
    bcftools: "\$(bcftools --version | awk '(NR==1){print \$NF}')"
EOF
"""
stub:
    def prefix = task.ext.prefix?:"${meta.id}.indep"
"""
touch versions.yml ${prefix}.bim ${prefix}.bed ${prefix}.fam 
"""
}
