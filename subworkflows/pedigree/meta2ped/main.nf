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
workflow META_TO_PED{
take:
	meta
	metas
main:
	MAKE_PED(
		meta,
		metas.map{[
			it.id,
			it.father?:"0",
			it.mother?:"0",
			it.sex?:"0",
			it.status?:"0"
			]}
		.map{it.join("\t")}
		.collect()
		)
	tmp = Channel.empty().ifEmpty([[id:"noped"],[]]).first()
emit:
	pedigree_gatk = tmp
	sample2collection = Channel.empty()
	versions = MAKE_PED.out.versions
}

process MAKE_PED {
input:
	val(meta)
	val(L)
output:
	tuple val(meta),path("*.gatk.ped"),optional:true,emit:gatk
	path("versions.yml"),emit:versions
script:
"""
cat << EOF | sort | uniq > raw.ped
${L.join("\n")}
EOF

python ${moduleDir}/ped.py raw.ped > /dev/null

touch versions.yml
"""

stub:


"""
touch versions.yml
"""
}
