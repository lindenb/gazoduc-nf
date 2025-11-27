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

process VISTA_DOWNLOAD {
tag "${meta.id?:dict.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../../conda/bioinfo.01.yml"
input:
    tuple val(meta),path(dict)
output:
	tuple val(meta),path("*.bed.gz"),path("*.bed.gz.tbi"),path("*.header"),emit:bed
	tuple val(meta),path("*.md"),emit:doc
	path("versions.yml"),emit:versions
script:
   	def TAG = "VISTA"
	def whatis="VISTA enhancers"
	def url =task.ext.url?:""
    if(url.isEmpty()) {
		if(meta.ucsc_name==null || meta.ucsc_name.isEmpty()) throw new IllegalArgumentException("${task.process} missing meta.ucsc name");
		url = "https://hgdownload.cse.ucsc.edu/gbdb/${meta.ucsc_name}/vistaEnhancers/vistaEnhancers.bb"
		}
	def jvm  = task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"
	def prefix = task.ext.prefix?:"${meta.id}.${TAG}"
"""
mkdir -p TMP/CACHE


if ${!url.isEmpty()}
then

	curl -L -o TMP/jeter.bb "${url}"

	bigBedToBed -udcDir=TMP/CACHE TMP/jeter.bb stdout |\\
			cut -f1,2,3,4 |\\
			jvarkit  ${jvm} bedrenamechr -R ${dict} --column 1 --convert SKIP  |\\
			LC_ALL=C sort  -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n -T TMP |\\
			bgzip > TMP/${prefix}.bed.gz




	echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="${whatis}">' > ${prefix}.header

else
	 
	touch  TMP/${prefix}.bed
	bgzip TMP/${prefix}.bed
	echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="Not available for ucsc: ${meta.ucsc_name}">' > ${prefix}.header

fi

tabix -p bed -f TMP/${prefix}.bed.gz
mv TMP/${prefix}.bed.gz ./
mv TMP/${prefix}.bed.gz.tbi ./

cat << EOF > ${prefix}.md
Vista enhancer.
EOF

cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${url}"
END_VERSIONS
"""

stub:
	def prefix = task.ext.prefix?:"vista"
"""
touch versions.yml ${prefix}.bed.gz ${prefix}.bed.gz.tbi ${prefix}.header ${prefix}.md
"""
}
