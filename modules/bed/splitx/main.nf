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

process BED_SPLITX {
tag "${meta.id} ${bed.name}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta),path(bed)
output:
    tuple val(meta),path("split*.bed"),optional:true,emit:beds
    path("versions.yml"),emit:versions
script:
    if(task.ext.njobs==null) throw new IllegalArgumentException("${task.process} task.n.jobs.undefined.")
    def njobs  = ((task.ext.njobs?:"10") as int)
    if(njobs<2) throw new IllegalArgumentException("${task.process} task.n.jobs<2")
"""
python3 ${moduleDir}/../../../src/python/split_bed.py" '${njobs}' 'split.${meta.id}' '${bed}' 

cat << EOF > versions.yml
${task.process}:
	python3: "\$(python3 --version | awk '{print \$2}')"
EOF
"""

stub:
    if(task.ext.njobs==null) throw new IllegalArgumentException("${task.process} task.n.jobs.undefined.")
"""
touch split.${meta.id}.0.bed split.${meta.id}.1.bed versions.yml
"""
}
