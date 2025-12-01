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

process ANNOTSV_APPLY {
tag "${meta.id?:""}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/annotsv.yml"
input:
    tuple val(meta1),path(annotsv_db)
    tuple val(meta2),path(dict) //for build/ucsc_name
    tuple val(meta),path(vcf)
output:
	tuple val(meta),path("*.annotSV.tsv"),optional:true,emit:annot
        tuple val(meta),path("*.unannotated.tsv"),optional:true,emit:unannotated
	path("versions.yml"),emit:versions
script:
   def prefix = task.ext.prefix?:"${meta.id}.annotsv"
   def ucsc_name = task.ext.ucsc_name?:"${meta2.ucsc_name}"
   def build = task.ext.build?:(ucsc_name=="hg38"?"GRCh38":(ucsc_name=="hg19"?"GRCh37":"unknown"))
   def args1 = task.ext.args1?:""
"""
mkdir -p TMP

# run annotSV
AnnotSV \\
	-annotationsDir "\${PWD}/${annotsv_db.name}" \\
	${args1} \\
	-genomeBuild ${build} \\
	-SVinputFile "${vcf}" \\
	-outputDir "TMP" \\
	-outputFile "${prefix}.annotSV.tsv"

mv -v TMP/*.tsv ./

touch versions.yml
"""

stub:
     def prefix = task.ext.prefix?:"${meta.id}.annotsv" 
"""
mkdir -p OUT
touch versions.yml OUT/${prefix}.annotSV.tsv OUT/${prefix}.annotSV.unannotated.tsv 
"""
}
