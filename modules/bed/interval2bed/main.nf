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

process INTERVAL_TO_BED {
executor "local"
tag "${contig}:${start1}-${end}"
input:
    tuple val(meta),val(contig),val(start1),val(end)
output:
    tuple val(meta),path("*.bed*"),emit:bed
script:
    def start0 = ((start1 as int) -1)
    def prefix = task.ext.prefix?:"${contig}_${start1}_${end}"
"""
echo  '${contig}\t${start0}\t${end}' >  ${prefix}.bed
"""

stub:
    def start0 = ((start1 as int) -1)
    def prefix = task.ext.prefix?:"${contig}_${start1}-${end}"
"""
echo  '${contig}\t${start0}\t${end}' >  ${prefix}.bed
"""
}
