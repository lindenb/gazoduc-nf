/*

Copyright (c) 2024 Pierre Lindenbaum

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


process CONCAT_FILES_01 {
afterScript "rm -rf TMP"
label "process_single"
input:
        path("INPUT/*")
output:
        path("concat*"),emit:output
script:
	def suffix = (task.ext && task.ext.suffix ? task.ext.suffix : ".list")
	def concat_n_files = (task.ext && task.ext.concat_n_files ? task.ext.concat_n_files:50)
	def downstream_cmd  = (task.ext && task.ext.downstream_cmd ? task.ext.downstream_cmd:"")
"""
hostname 1>&2

mkdir -p TMP
find ./INPUT/ -type l > TMP/jeter.txt

xargs -a TMP/jeter.txt -L ${concat_n_files} cat ${downstream_cmd} > TMP/concat.data

mv TMP/concat.data "concat${suffix}"

sleep 5
"""
}
