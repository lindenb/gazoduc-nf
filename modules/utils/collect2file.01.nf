/*

Copyright (c) 2022 Pierre Lindenbaum

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

process COLLECT_TO_FILE_01 {
executor "local"
tag "N=${L.size()}"
input:
	val(meta)
	val(L)
output:
	path("concat.list"),emit:output
script:
"""
hostname 1>&2
set -o pipefail

cat << EOF | awk -F '/' '{printf("%s\t%s\\n",\$NF,\$0);}' | sort -t '\t' -T. -k1,1 -k2,2 | cut -f 2 | uniq > concat.list
${L.join("\n")}
EOF

## if it's too fast, prevent clock problems
sleep 10
touch -c concat.list
"""
stub:
"""
touch concat.list
"""
}
