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
include {isBlank      } from  '../../../modules/utils/functions.nf'

process DOWNLOAD_GNOMAD_SV {
tag "${meta1.id?:fasta.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
output:
	tuple val(meta1),path("*.vcf.gz"),path("*.tbi"),emit:vcf
	path("versions.yml"),emit:versions
script:
	def url = task.ext.url?:"";
    def prefix = task.ext.suffix?:""
    if(!isBlank(url)) {
        //nothing
        }
    else if(isBlank(meta1.ucsc_name)) {
        url = ""
        }
    else if(meta1.ucsc_name == "hg19") {
        url = "https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz";
        if(isBlank(prefix)) prefix = "gnomad_v2.1_sv.sites"
        }
    else if(meta1.ucsc_name == "hg38") {
        url = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz";
        if(isBlank(prefix)) prefix = "gnomad_v4.1_sv.sites"
        }
    else {
        url = ""
        }
    if(isBlank(prefix)) prefix="gnomadsv.sites";
	def suffix =  task.ext.suffix?:".vcf.gz"
    def jvm = task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData"
    if(isBlank(url)) {
        log.warn("Empty url for ${task.process}");
        }
"""
hostname 1>&2
mkdir -p TMP

if ${!isBlank(url)}
then
    curl -L -o TMP/jeter.vcf.gz -L "${url}"
else

# create dummy VCF, in the telomeric region of chr22

cat << EOF > jeter.vcf
#fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LOW_CALL_RATE,Description="">
##INFO=<ID=POPMAX_AF,Number=1,Type=Float,Description="Maximum allele frequency across any population (biallelic sites only).">
22\t1\t.\tN\t<DEL>\t1\tLOW_CALL_RATE \tPOPMAX_AF=0.9
EOF

bgzip jeter.vcf

fi

bcftools query -f "%CHROM\\n"  TMP/jeter.vcf.gz |\\
    uniq | sort -T TMP | uniq |\\
    awk '{printf("%s\t%s\\n",\$1,\$1);}' |\\
	jvarkit  bedrenamechr --column 2 -R "${fasta}" --convert SKIP > TMP/chroms.tsv

	bcftools annotate --rename-chrs TMP/chroms.tsv TMP/jeter.vcf.gz |\\
	bcftools sort \\
        -T TMP/sort \\
        --max-mem "${task.memory.giga}G" \\
        -O z \\
        -o TMP/jeter.vcf.gz

    bcftools index \\
        --threads ${task.cpus} \\
        --force \\
        --tbi \\
        TMP/jeter.vcf.gz

mv  TMP/jeter.vcf.gz "${prefix}.vcf.gz"
mv  TMP/jeter.vcf.gz.tbi "${prefix}.vcf.gz.tbi"

cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${url}"
END_VERSIONS
"""

stub:
	def prefix = "${meta1.id}"
"""
touch versions.yml "${prefix}.vcf.gz" "${prefix}.vcf.gz.tbi"
"""
}
