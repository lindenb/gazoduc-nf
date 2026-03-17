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
include { parseBoolean  } from '../functions.nf'
include { verify        } from '../functions.nf'
include { isBlank       } from '../functions.nf'

process SIMPLE_COMMAND {
tag "${meta.id}"
label "process_single"
executor "local"
afterScript "rm -rf TMP"
input:
    tuple val(meta),path(input)
output:
    tuple val(meta),path("${meta.prefix?:""}concat${meta.suffix?:".txt"}"),emit:output
	path("versions.yml"),emit:versions
script:
	def command = task.ext.command?:""
	verify(!isBlank(command),"${task.process} task.ext.command is blank")
	def downstream = meta.downstream?:"| cat "
	def prefix = task.ext.prefix?:"${meta.id}"
"""
hostname 1>&2

${command} ${input} ${downstream} > '${meta.prefix?:""}concat${meta.suffix?:".txt"}'


# avoid timestamp problem
sleep 2

touch versions.yml
"""
stub:
"""
touch versions.yml
"""
}

