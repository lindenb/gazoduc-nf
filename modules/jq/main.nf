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

process JQ {
tag "${meta.id?:json.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../conda/jq.yml"
input:
    tuple val(meta ),path(json)
    tuple val(meta2),path(opt_query)
output:
	tuple val(meta),path("*.txt"),emit:output
	path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"${meta.id}"
    def query0 = task.ext.query?:""
    def query_str  =  (opt_query?"":query0)
    def args1  = task.ext.args1?:""
    def post_cmd = task.ext.post_cmd?:""
    if((opt_query?false:true) && query_str.trim().isEmpty()) throw new IllegalArgumentException("${task.process} query undefined");
"""
mkdir -p TMP

if ${opt_query?true:false}
then
    jq ${args1} --from-file "${opt_query}" "${json}" ${post_cmd} > TMP/jeter.out.txt
else
    jq ${args1} '${query_str}' "${json}" ${post_cmd} > TMP/jeter.out.txt
fi

mv TMP/jeter.out.txt ${prefix}.out.txt

cat << END_VERSIONS > versions.yml
"${task.process}":
	jq: "todo"
END_VERSIONS
"""
}
