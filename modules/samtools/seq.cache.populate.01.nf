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

include {getVersionCmd;moduleLoad} from '../utils/functions.nf'

/** apply seq_cache_populate.pl to a reference . see http://www.htslib.org/doc/samtools.html 

then use:

	export REF_CACHE=${refcache}/%2s/%2s/%s

*/
process SEQ_CACHE_POPULATE_01 {
tag "${reference}"
afterScript "rm -rf TMP"
input:
      	val(meta)
        val(reference)
output:
	path("CACHE"),emit:cache
	path("version.xml"),emit:version
script:
	def subdirs = meta.seq_cache_populate_subdirs?:"2"
"""
hostname 1>&2
${moduleLoad("samtools")}
mkdir -p TMP

seq_cache_populate.pl -root \${PWD}/TMP  -subdirs "${subdirs}" "${reference}"

mv TMP CACHE

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">invoke samtools/seq_cache_populate.pl to build REF_CACHE. See <a>http://www.htslib.org/doc/samtools.html</a>.</entry>
        <entry key="reference">${reference}</entry>
        <entry key="subdirs">${subdirs}</entry>
        <entry key="versions">${getVersionCmd("samtools")}</entry>
</properties>
EOF
"""
}

