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

process HMC_DOWNLOAD{
tag "${meta1.id?:fasta.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
output:
	tuple val(meta1),path("*.bed.gz"), path("*.bed.gz.tbi"), path("*.header"),emit:bed
	path("versions.yml"),emit:versions
script:
	def TAG=  task.ext.tag?:"HMC"
	def url="";
    if(meta1.ucsc_name && (meta1.ucsc_name.equals("hg19") || meta1.ucsc_name.equals("hg38"))) {
        url = "http://hgdownload.soe.ucsc.edu/gbdb/${meta1.ucsc_name}/hmc/hmc.bw"
        }
    if(url.trim().isEmpty()) {
        throw new IllegalArgumentException("bad url for ${task.process}");
    }
"""
hostname 1>&2
mkdir -p TMP

curl -L -o TMP/jeter.bw "${url}"

bigWigToBedGraph TMP/jeter.bw stdout |\
	jvarkit   -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP bedrenamechr -f "${fasta}" --column 1 --convert SKIP  |\
		LC_ALL=C sort  -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n |\
		bgzip > TMP/${TAG}.bed.gz

tabix -p bed -f TMP/${TAG}.bed.gz


mv TMP/${TAG}.bed.gz ./
mv TMP/${TAG}.bed.gz.tbi ./

echo '##INFO=<ID=${TAG},Number=1,Type=Float,Description="Homologous Missense Constraint (HMC) is a amino acid level measure of genetic intolerance of missense variants within human populations. For all assessable amino-acid positions in Pfam domains, the number of missense substitutions directly observed in gnomAD (Observed) was counted and compared to the expected value under a neutral evolution model (Expected). The upper limit of a 95% confidence interval for the Observed/Expected ratio is defined as the HMC score. Missense variants disrupting the amino-acid positions with HMC<0.8 are predicted to be likely deleterious. This score only covers PFAM domains within coding regions . URL= ${url}">' > ${TAG}.header

cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "\${URL}"
END_VERSIONS
"""
}
