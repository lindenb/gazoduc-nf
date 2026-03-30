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

/*

Count the number of variants in each chunk
Column 3: total number of variants
Column 4: comma-deliminted number of sub-chunks after accounting for size
Column 5: comma-deliminted number of sub-chunks after accounting for size and distance again

positional arguments:
  input                 Input VCF

options:
  -h, --help            show this help message and exit
  -b BED, --bed BED     Bed file of variants to chunk
  -o OUTPUT, --output OUTPUT
                        Output name
  -c CHUNKSIZE, --chunksize CHUNKSIZE
                        Distance between variants to split chunks (500)
  -s SIZEMIN, --sizemin SIZEMIN
                        Minimum SV length (50)
  -S SIZEMAX, --sizemax SIZEMAX
                        Maximum SV length (50000)


*/

process TRUVARI_ANNO_CHUNCKS {
    label "process_single"
	  tag "${meta.id}"
    afterScript "rm -rf TMP"
    conda "${moduleDir}/../../../conda/truvari.01.yml"
    input:
      tuple val(meta1),path(opt_bed) //Bed file of variants to chunk
		  tuple val(meta ),path(vcf)//generated with bcftools merge
   output:
		tuple val(meta),path("*.bed"),emit: bed
		path("versions.yml"),emit:versions
  script:
		def args1 = task.ext.args1?:""
    def size_min = task.ext.size_min?:50
    def size_max = task.ext.size_max?:1000000
    def prefix= task.ext.prefix?:"${meta.id}.anno_chunks"
    """
	hostname 1>&2
	mkdir -p TMP
	
     truvari anno chunks \\
	${args1} \\
        ${opt_bed?"--bed \"${opt_bed}\" ":""} \\
        --sizemin ${size_min} \\
        --sizemax ${size_max} \\
	${vcf} > ${prefix}.bed

cat << END_VERSIONS > versions.yml
"${task.process}":
	truvari: \$(truvari version)
END_VERSIONS
"""

stub:
def prefix= task.ext.prefix?:"${meta.id}.anno_chunks"
"""
touch versions.yml ${prefix}.bed
"""
}
