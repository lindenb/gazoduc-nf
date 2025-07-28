

process BCFTOOLS_CALL {
tag "${meta.id?:""} ${optional_bed?optional_bed.name:""}"
label "process_single"
afterScript "rm -rf TMP TMP2"
array 100
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(optional_pedigree)
    tuple val(meta4),path(optional_ploidy)
    tuple val(meta),path("BAMS/*"),path(optional_bed)
output:
    tuple val(meta),path("*.bcf",arity:'1'),path("*.csi",arity:'1'),path(optional_bed),emit:vcf
output
    path("versions.yml"),emit:versions
	def ploidy = "";
      if(optional_ploidy?true:false) {
        ploidy = "--ploidy ${optional_ploidy.name}";
        }
    else if(meta1.ucsc_name && meta1.ucsc_name.equals("hg19")) {
            ploidy = "--ploidy GRCh37";
    } else if(meta1.ucsc_name && meta1.ucsc_name.equals("hg19")) {
        ploidy = "--ploidy GRCh38";
    } 

script:
    def args1 = task.ext.args1?:""
    def args2 = task.ext.args2?:""
    def prefix = task.ext.prefix?:"\${MD5}"+(optional_bed?"."+optional_bed.baseName:"")
    def mapq = task.ext.mapq?:10
"""
set -x
mkdir -p TMP
find BAMS/ -name "*am" | samtools samples | sort -T TMP -k1,1 -k2,2 -t '\t' | cut -f 2 > TMP/bams.list
test -s  TMP/bams.list
MD5=`cat TMP/bams.list | md5sum | cut -d ' ' -f1`


bcftools mpileup --redo-BAQ -a 'FORMAT/AD' -a 'FORMAT/DP' -a 'INFO/AD' \\
	--fasta-ref "${fasta}" \\
	--threads ${task.cpus} \\
	${optional_bed?"--regions-file \"${optional_bed}\"":""}" \\
	-q ${mapq} \\
	-O u \\
    ${args1} \\
	--bam-list TMP/bams.list \\
	-o TMP/jeter.bcf

bcftools call  \\
	--keep-alts \\
	${ploidy} \\
	--threads ${task.cpus} \\
	--multiallelic-caller \\
	--variants-only \\
	--output-type b \\
	-o TMP/jeter2.bcf \\
    ${pedigree?"--samples-file ${pedigree}":""} \\
    ${args2} \\
	TMP/jeter.bcf

bcftools  +fill-tags --threads ${task.cpus}  -O b  -o TMP/jeter.bcf  TMP/jeter2.bcf -- -t AN,AC,AF,AC_Hom,AC_Het,AC_Hemi,NS
bcftools index  --threads ${task.cpus}  TMP/jeter.bcf

mv TMP/jeter.bcf     ${prefix}.bcf
mv TMP/jeter.bcf.csi ${prefix}.bcf.csi

cat << END_VERSIONS > versions.yml
${task.process}:
    bcftools: \$(bcftools version | awk '(NR==1)  {print \$NF}')
END_VERSIONS
"""
}
