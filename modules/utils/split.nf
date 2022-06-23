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
/** use linux split to split a file into parts */
process LINUX_SPLIT {
executor "local"
tag "${file(filein).name}"
input:
      	val(meta)
        val(filein)
output:
        path("clusters.list"),emit:output
	path("version.xml"),emit:version
script:
       	def prefix = meta.prefix?:"chunck."
       	def suffix = meta.suffix?:".txt"
       	def method = meta.method?:"--lines=10"
        """

	mkdir TMP
        split -a 9 --additional-suffix=${suffix} ${method} "${filein}" "TMP/${prefix}"

      	find \${PWD}/TMP -type f -name "${prefix}*${suffix}" > clusters.list

	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="Name">${task.process}</entry>
		<entry key="Description">Split file into parts</entry>
		<entry key="Input">${filein}</entry>
		<entry key="method">${method}</entry>
		<entry key="N-FILES">\$(wc -l < clusters.list)</entry>
		<entry key="Output">clusters.list</entry>
	</properties>
	EOF
        """
	}
