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
process PER_CONTIG {
tag "${meta.id}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
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
    tuple val(meta ),path("indep_*"),emit:bfile
    path("versions.yml"),emit:versions
script:
    def norm_contig = meta.contig
    def optional_contig_1k = metaK.contig
    def optional_contig_gnomad = metaG.contig
    def arg1 = ""
    def cpus3 = task.cpus>3?"--thread ${(task.cpus/3) as int}":""
    def plink_args  = "--const-fid 1 --allow-extra-chr --allow-no-sex --threads ${task.cpus}"
    def prefix = task.ext.prefix?:"${meta.id}"
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

bcftools view -G --no-header TMP/jeter1.bcf |head 1>&2 || true
bcftools query -l TMP/jeter1.bcf | cat -n | tail 1>&2 || true


# 
# 1000 GENOMES DATA
#
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
			jvarkit -Djava.io.tmpdir=TMP bedrenamechr -f "${gnomad}" --column 1 --convert SKIP |\\
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
			jvarkit  -Djava.io.tmpdir=TMP bedrenamechr -f "${fasta}" --column 1 --convert SKIP |\\
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

	bcftools query -f '.\\n' TMP/jeter1.bcf | wc -l 1>&2
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

bcftools query -f '.\\n' TMP/jeter1.bcf | wc -l 1>&2

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

bcftools query -f '.\\n' TMP/jeter1.vcf.gz | wc -l 1>&2
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
    plink: "\$(plink --version | awk '{print \$2}')"
    bcftools: "\$(bcftools --version | awk '(NR==1){print \$NF}')"
EOF
"""
stub:
def prefix = "indep_${norm_contig}"
"""
touch versions.yml ${prefix}.bim ${prefix}.bed ${prefix}.fam 
"""
}
