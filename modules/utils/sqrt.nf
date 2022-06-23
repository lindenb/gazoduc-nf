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

process SQRT_FILE {
executor "local"
tag "${file(files).name}"
input:
      	val(meta)
        val(files)
output:
        path("clusters.list"),emit:clusters
	path("version.xml"),emit:version
script:
       	def min_file_split = getKeyValue(meta,"min_file_split","-1")
        """
        SQRT=`awk 'END{X=NR;if(${min_file_split} > 0 && X <= ${min_file_split}){print(X);} else {z=sqrt(X); print (z==int(z)?z:int(z)+1);}}' "${files}"`

        split -a 9 --additional-suffix=.list --lines=\${SQRT} "${files}" chunck.

      	find \${PWD} -type f -name "chunck.*.list" > clusters.list

	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="Name">${task.process}</entry>
		<entry key="Description">Split file into parts</entry>
		<entry key="Input">${files}</entry>
		<entry key="N">\${SQRT}</entry>
		<entry key="N-FILES">\$(wc -l < clusters.list)</entry>
		<entry key="Output">clusters.list</entry>
	</properties>
	EOF
        """
	}
