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
process XHUNTER_DENOVO_CASE_CONTROL {
    label "process_single"
	tag "${meta.id}"
	afterScript "rm -rf TMP"
    conda "${moduleDir}/../../../conda/multiqc.yml"//yes not xhunter, just need numpy
	input:
		tuple val(meta1),path(fasta)
        tuple val(meta2),path(fai)
		tuple val(meta3),path(json)
		tuple val(meta ),path(manifest)
	output:
		tuple val(meta),path("*.tsv"),emit:tsv
		path("versions.yml"),emit:versions
	script:
        def method = task.ext.method?:""
        if(method.trim().isEmpty()) throw new IllegalArgumentException("XHUNTER_DENOVO_CASE_CONTROL task.ext.method is empty");
        def prefix=task.ext.prefix?:"${meta.id}.${method}"
	"""
	hostname 1>&2
	mkdir -p TMP
	curl -L -o TMP/jeter.zip "https://github.com/Illumina/ExpansionHunterDenovo/archive/refs/heads/master.zip"
	(cd TMP && unzip jeter.zip)

	TMP/ExpansionHunterDenovo-master/scripts/casecontrol.py  ${method} \
		--manifest '${manifest}' \
		--multisample-profile '${json}' \
		--output TMP/jeter.tsv 


	mv TMP/jeter.tsv "${prefix}.tsv"

cat << EOF > versions.yml
"${task.process}":
    xhunter: todo
EOF
	"""


stub:
        
def method = task.ext.method?:""
if(method.trim().isEmpty()) throw new IllegalArgumentException("XHUNTER_DENOVO_CASE_CONTROL task.ext.method is empty");
"""
touch versions.yml ${meta.id}.${method}.tsv
"""
}

