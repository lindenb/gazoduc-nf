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
process BCFTOOLS_CALL {
tag "${meta.id?:""} ${optional_bed?optional_bed.name:""}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(optional_pedigree)
    tuple val(meta4),path(optional_ploidy)
    tuple val(meta ),path("BAMS/*"),path(optional_bed)
output:
    tuple val(meta),path("*.bcf",arity:'1'),path("*.csi",arity:'1'),path(optional_bed),emit:vcf
    path("versions.yml"),emit:versions
script:
    def args1 = task.ext.args1?:""
    def args2 = task.ext.args2?:""
    def prefix = task.ext.prefix?:(meta.id?meta.id+".":"")+"\${MD5}"+(optional_bed?"."+optional_bed.baseName:"")
    def mapq = task.ext.mapq?:10

    def ploidy = "";
    if(optional_ploidy?true:false) {
        ploidy = "--ploidy ${optional_ploidy.name}";
        }
    else if(meta1.ucsc_name && meta1.ucsc_name.equals("hg19")) {
            ploidy = "--ploidy GRCh37";
    } else if(meta1.ucsc_name && meta1.ucsc_name.equals("hg38")) {
        ploidy = "--ploidy GRCh38";
    } 
"""
set -x
mkdir -p TMP
find BAMS/ -name "*am" | samtools samples | sort -T TMP -k1,1 -k2,2 -t '\t' > TMP/jeter.samples.tsv

cut -f 2 TMP/jeter.samples.tsv > TMP/bams.list
test -s  TMP/bams.list
MD5=`cat TMP/bams.list | md5sum | cut -d ' ' -f1`

if ${optional_pedigree?true:false}
then

    tr " " "\t" < "${optional_pedigree}" | tr -s "\t" |  sort -T TMP -k2,2 -t '\t' > TMP/jeter2.tsv

    join -t '\t' -1 1 -2 2 -o '2.1,2.2,2.3,2.4,2.5' TMP/jeter.samples.tsv TMP/jeter2.tsv |\\
        awk -F '\t' '{OFS="\t";if(\$5=="male" || \$5=="XY"){\$5="1";}else if(\$5=="female" || \$5=="XX"){\$5="2";} else {\$5="1";} print;}' > TMP/jeter.pedigree

fi

bcftools mpileup --redo-BAQ -a 'FORMAT/AD' -a 'FORMAT/DP' -a 'INFO/AD' \\
	--fasta-ref "${fasta}" \\
	--threads ${task.cpus} \\
	${optional_bed?"--regions-file \"${optional_bed}\"":""} \\
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
	--output-type u \\
	-o TMP/jeter2.bcf \\
    ${optional_pedigree?"--samples-file TMP/jeter.pedigree":""} \\
    ${args2} \\
	TMP/jeter.bcf

bcftools sort \\
    --max-mem '${task.memory.giga}G' \\
    -T TMP/sort \\
    -O u \\
    -o TMP/jeter.bcf \\
    TMP/jeter2.bcf

mv  TMP/jeter.bcf TMP/jeter2.bcf

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
