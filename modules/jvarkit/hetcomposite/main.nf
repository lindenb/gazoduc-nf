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

process HET_COMPOSITE {
tag "${meta.id?:""}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta4),path(pedigree)
    tuple val(meta),path(vcf),path(vcfidx) // MUST BE ANNOTATED WITH SNPEFF, REMOVE FREQUENT VARIANTS
output:
	tuple val(meta),path("*.bcf"),path("*.bcf.csi"),emit:vcf
	tuple val(meta),path("*.genes.report"),emit:genes_report
	tuple val(meta),path("*.variants.report"),emit:variants_report
	path("versions.yml"),emit:versions
when:
	task.ext.when == null || task.ext.when
script:
    def prefix = task.ext.prefix?:vcf.baseName+".hetcomposite"
	def extractors = task.ext.extractors?:"ANN/GeneId"
	def max_variants = task.ext.max_variants?:30
"""
hostname 1>&2

mkdir -p TMP
set -x

# convert pedigree if no 6th column
awk '{S=\$6 ; if(NF==5 || S=="") { if(\$3!="0" && \$4!="0") {S="case";} else {S="control"} }  printf("%s\t%s\t%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,\$4,\$5,S);}' ${pedigree} > TMP/pedigree.tsv

## all other samples are controls
comm -13 \\
	<(cut -f 2 TMP/pedigree.tsv  | sort | uniq) \\
	<(bcftools query -l '${vcf}'| sort | uniq) |\\
	awk '{printf("if(!acceptControl(variant,\\"%s\\")) return false;\\n",\$1);}' >> TMP/custom.m4

awk -F '\t' '(\$6=="control" || \$6=="unaffected") {printf("if(!acceptControl(variant,\\"%s\\")) return false;\\n",\$2);}'  TMP/pedigree.tsv >> TMP/custom.m4

awk -F '\t' '(\$6=="case" || \$6=="affected") {printf("if(acceptTrio(variant,\\"%s\\",\\"%s\\",\\"%s\\")) return true;\\n",\$2,\$3,\$4);}' TMP/pedigree.tsv  >> TMP/custom.m4


m4 -P -I TMP < "${moduleDir}/select.m4"  > TMP/jeter.code

bcftools view "${vcf}"  |\\
	jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP vcffilterjdk --body -f TMP/jeter.code > TMP/jeter2.vcf
mv TMP/jeter2.vcf TMP/jeter1.vcf


jvarkit -Xmx${task.memory.giga}G  -Djava.io.tmpdir=TMP vcfcomposite \\
	--extractors "${extractors}" \\
	--filter "" \\
	--genes ${prefix}.genes.report \\
	--pedigree  TMP/pedigree.tsv \\
	--report ${prefix}.variants.report \\
	--tmpDir TMP \\
	--max-variants ${max_variants} \\
	TMP/jeter1.vcf > TMP/jeter2.vcf

mv TMP/jeter2.vcf TMP/jeter1.vcf


bcftools sort -T TMP/sort -o TMP/${prefix}.bcf -O b TMP/jeter1.vcf
bcftools index --threads ${task.cpus} TMP/${prefix}.bcf

mv TMP/${prefix}.bcf ./
mv TMP/${prefix}.bcf.csi ./

cat << EOF > versions.yml
${task.process}:
	jvarkit: TODO
	max-variants: "${max_variants}"
	extractors: "${extractors}"
EOF
"""
}

