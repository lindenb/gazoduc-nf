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
process GET_MAPPABILITY {
label "process_single"
tag "${fasta.name}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(fai)
output:
	tuple path("blacklist.*"),optional:true,emit:output
	path("versions.yml"),emit:versions
script:
	def url0 = meta.ucsc_name.equals("hg19")
			? "https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh37.dna.primary_assembly.fa.r101.s501.blacklist.gz"
			: (meta.ucsc_name.equals("hg38")
			? "https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz"
			: "")
	def url = task.ext.url?:url0
	def prefix = task.ext.prefix?:"blacklist"
"""
mkdir -p TMP


if test ! -z "\${url}"
then
	wget -O TMP/${prefix}.gz "\${URL}"
	wget -O TMP/${prefix}.gz.fai "\${URL}.fai"
	wget -O TMP/${prefix}.gz.gzi "\${URL}.gzi"
	mv TMP  mappability
fi

"""
}
