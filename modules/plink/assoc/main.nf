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

process PLINK_ASSOC {
tag "${meta.id}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1 ),path(mds)
	tuple val(meta ),path(bim),path(bed),path(fam)
output:
    tuple val(meta ),path("*.qassoc", arity: '1..*'),emit:assoc
	tuple val(meta ),path("*.log"),emit:log
    path("versions.yml"),emit:versions
script:
    def plink_args  = "--const-fid 1 --allow-extra-chr --allow-no-sex "
	def prefix = task.ext.prefix?:"${meta.id}"
	"""
	mkdir -p TMP

	# create pheno file https://www.cog-genomics.org/plink/1.9/input#pheno
	awk '(NR==1) {split(\$0,header);} {X=0;for(i=1;i<=NF;i++) {if(header[i] !="SOL") {printf("%s%s",(X==0?"":"\t"),\$i);X=1;}} printf("\\n");}' '${mds}' > TMP/pheno.tsv

	plink \\
        ${plink_args} \\
		--threads ${task.cpus} \\
        --bfile ${bim.baseName} \\
		--pheno TMP/pheno.tsv \\
		--all-pheno \\
		--assoc \\
		--out TMP/${prefix}
	
mv TMP/${prefix}* ./

cat << EOF > versions.yml
${task.process}:
    plink: \$(plink --version | awk '(NR==1) {print \$2}')
EOF
"""
stub:
"""
touch versions.yml
"""
}
