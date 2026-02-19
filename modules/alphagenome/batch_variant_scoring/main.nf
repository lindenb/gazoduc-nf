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
include {verify   } from '../../../modules/utils/functions'

process ALPHAGENOME_BATCH_VARIANT_SCORING {
tag "${meta.id?:""}"
afterScript "rm -rf TMP"
label "process_single"
secret 'ALPHAGENOME_API_KEY'
maxForks 1
conda "${moduleDir}/../../../conda/alphagenome.yml"
input:
    tuple val(meta ),path(variants) /* tsv file variant_id\tCHROM\POS\REF\ALT */
output:
    tuple val(meta),path("*.tsv.gz"),emit:tsv
script:
    def organism  = task.ext.organism?:"human"
    def sequence_length = task.ext.sequence_length ?:"1MB" 
	def prefix = task.ext.prefix?:"${meta.id}.variant_scores"
    // https://www.alphagenomedocs.com/api/models.html#variant-scorers
    //def scorers = task.ext.scorers?:"variant_scorers.RECOMMENDED_VARIANT_SCORERS"
"""
hostname 1>&2

ulimit -c unlimited

python3 "${moduleDir}/script.py" \\
    "\${ALPHAGENOME_API_KEY}" \\
    "${variants}" \\
    '${organism}' \\
    '${sequence_length}'  

gzip --best variant_scores.tsv
mv variant_scores.tsv.gz ${prefix}.tsv.gz

cat << END_VERSIONS > versions.yml
"${task.process}":
	alphagenome: ""
END_VERSIONS
"""
stub:
	def prefix = task.ext.prefix?:"${meta.id}.variant_scores" 
"""
touch versions.yml ${prefix}.tsv
gzip ${prefix}.tsv
"""
}
