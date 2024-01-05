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
include { getKeyValue; getClassTaxonomy } from './functions.nf'

/** break a file into set of SQRT(count-lines(file)) files */
process SQRT_FILE {
tag "${filein.name}"
executor "local"
input:
      	val(meta)
        path(filein)
output:
        path("clusters${meta.suffix?:".list"}"),emit:output
	path("version.xml"),emit:version
script:
	if(!meta.containsKey("suffix")) throw new IllegalArgumentException("suffix missing");
	if(!meta.containsKey("min_file_split")) throw new IllegalArgumentException("min_file_split missing");
	def suffix = meta.suffix?:".list"
       	def min_file_split = meta.min_file_split?:-1
        """
        SQRT=`awk 'END{X=NR;if(${min_file_split} > 0 && X <= ${min_file_split}){print(X);} else {z=sqrt(X); print (z==int(z)?z:int(z)+1);}}' "${filein}"`
	mkdir -p OUT
        split -a 9 --additional-suffix=${suffix} --lines=\${SQRT} "${filein}" OUT/chunck.

      	find \${PWD}/OUT -type f -name "chunck.*${suffix}" > clusters.list

	# if too fast, prevent problems with timestamp ?
	sleep 10
	touch -c clusters.list

	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="Name">${task.process}</entry>
		<entry key="Description">Split file into parts</entry>
		<entry key="Input">${filein}</entry>
		<entry key="N">\${SQRT}</entry>
		<entry key="N-FILES">\$(wc -l < clusters.list)</entry>
		<entry key="Output">clusters.list</entry>
	</properties>
	EOF
        """
	}

