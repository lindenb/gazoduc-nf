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
	tuple val(meta1),path("*.{vcf.gz,bcf}"),path("*.{csi,tbi}"),emit:vcf
	path("versions.yml"),emit:versions
script:
	def url = "";
    if(meta1.ucsc_name && meta1.ucsc_name.equals("hg19")) {
        url = "https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz";
        }
    else if(meta1.ucsc_name && meta1.ucsc_name.equals("hg38")) {
        url = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz";
        }
    else {
        throw new IllegalArgumentException("${task.process}: url is undefined.");
        }
    def prefix  = task.ext.prefix?:"gnomadsv.sites"
	def suffix =  task.ext.suffix?:".bcf"
"""
hostname 1>&2
mkdir -p TMP

curl -o TMP/jeter.vcf.gz -L "${url}"

bcftools query -f "%CHROM\\n"  TMP/jeter.vcf.gz |\\
    uniq | sort -T TMP | uniq |\\
    awk '{printf("%s\t%s\\n",\$1,\$1);}' |\\
	jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData bedrenamechr --column 2 -R "${fasta}" --convert SKIP > TMP/chroms.tsv

	bcftools annotate --rename-chrs TMP/chroms.tsv TMP/jeter.vcf.gz |\\
	bcftools sort \\
        -T TMP/sort \\
        --max-mem "${task.memory.giga}G" \\
        -O ${suffix.contains("b")?"b":"z"} \\
        -o TMP/jeter.bcf

    bcftools index \\
        --threads ${task.cpus} \\
        --force \\
        ${suffix.contains("b")?"":"--tbi"} \\
        TMP/jeter.bcf

mv  TMP/jeter.bcf "${prefix}${suffix}"
mv  TMP/jeter.bcf.csi "${prefix}${suffix}.csi"

cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${url}"
END_VERSIONS
"""
}
