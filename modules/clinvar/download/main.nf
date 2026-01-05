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


process CLINVAR_DOWNLOAD {
afterScript "rm -rf TMP"
tag "${meta.id}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta),path(dict)
output:
    tuple val(meta),path("*.bcf"),path("*.bcf.csi"),optional:true,emit:vcf
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"${meta.id}.clinvar"
    def base = task.ext.url?:"https://ftp.ncbi.nlm.nih.gov/pub/clinvar"
    def url = task.ext.url?:""
    def ucsc_name = task.ext.ucsc_name?:meta.ucsc_name
    if(!url.isEmpty()) {
	// nothing
	}
    else if(ucsc_name!=null && ucsc_name=="hg38") {
	url = "${base}/vcf_GRCh38/clinvar.vcf.gz"
	}
    else if(ucsc_name!=null && ucsc_name=="hg19") {
	url = "${base}/vcf_GRCh37/clinvar.vcf.gz"
	}
    def local_vcf = task.ext.local_vcf?:"NO_FILE" 

"""
set -o pipefail
set -x
mkdir -p TMP


if test -f "${local_vcf}" && test -f "${local_vcf}.tbi"
then
 
bcftools view -O v "${local_vcf}" |\\
        jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  vcfsetdict -R "${dict}"  -n SKIP |\\
        bcftools view -O b -o TMP/jeter.bcf

elif ${!url.isEmpty()}
then

curl -L  "${url}" |\\
        bcftools view -O v  |\\
        jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  vcfsetdict -R "${dict}"  -n SKIP |\\
        bcftools view -O b -o TMP/jeter.bcf

fi


if test -f TMP/jeter.bcf
then

bcftools view --header-only TMP/jeter.bcf | grep "^##INFO" | cut -d '=' -f 3 | cut -d, -f 1| grep -v "#"  |\\
        awk '{printf("INFO/%s\tCLINVAR_%s\\n",\$1,\$1);}' > TMP/rename.tsv

bcftools annotate \\
    --threads ${task.cpus} \\
    --rename-annots TMP/rename.tsv \\
    -O b \\
    -o TMP/${prefix}.bcf \\
    TMP/jeter.bcf

bcftools index  \\
    --threads ${task.cpus} \\
     TMP/${prefix}.bcf

mv TMP/${prefix}.bcf ./
mv TMP/${prefix}.bcf.csi ./

fi


cat << EOF > versions.yml
${task.process}:
    bcftools: TODO
EOF
"""

stub:
    def prefix = task.ext.prefix?:"${meta.id}.clinvar"
"""
touch versions.yml ${prefix}.bcf ${prefix}.bcf.csi
"""
}
