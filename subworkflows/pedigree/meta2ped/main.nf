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
		metas
			.map{[ /* weird bug, if I don't use quote, just it?:something, it doesn't work... */
			"${it.id?:(it.sample?:"")}",
			"${it.father?:""}",
			"${it.mother?:""}",
			"${it.sex?:""}",
			"${it.status?:""}",
			"${it.collection?:(it.population?:"")}"
			]}
		.map{it.join("\t")}
		.collect()
		)

	cases = MAKE_PED.out.cases.ifEmpty([[id:"no.cases"],[]])
	controls = MAKE_PED.out.cases.ifEmpty([[id:"no.controls"],[]])
	pedigree = MAKE_PED.out.ped.ifEmpty([[id:"noped"],[]])
	sample2collection = MAKE_PED.out.sample2col.ifEmpty( [[id:"nosn2col"],[]])
	pedigree_gatk  = MAKE_PED.out.gatk.ifEmpty( [[id:"nogatkped"],[]])
	sample2status = MAKE_PED.out.sample2status.ifEmpty( [[id:"nosample2status"],[]])
emit:
	cases
	controls
	pedigree
	pedigree_gatk
	sample2collection
	sample2status
	versions = MAKE_PED.out.versions
}

process MAKE_PED {
label "process_single"
conda "${moduleDir}/../../../conda/multiqc.yml"
input:
	val(meta)
	val(L)
output:
	tuple val(meta),path("pedigree4gatk.ped"),optional:true,emit:gatk
	tuple val(meta),path("raw.ped"),optional:true,emit:ped
	tuple val(meta),path("males.txt"),optional:true,emit:males
	tuple val(meta),path("females.txt"),optional:true,emit:females
	tuple val(meta),path("cases.txt"),optional:true,emit:cases
	tuple val(meta),path("controls.txt"),optional:true,emit:controls
	tuple val(meta),path("sample2collection.tsv"),optional:true,emit:sample2col
	tuple val(meta),path("sample2status.tsv"),optional:true,emit:sample2status
	
	path("versions.yml"),emit:versions
script:
"""
cat << EOF | sort -T . | uniq > raw.ped
${L.join("\n")}
EOF

python3 ${moduleDir}/ped.py raw.ped > /dev/null

for F in males.txt females.txt cases.txt controls.txt sample2collection.tsv pedigree4gatk.ped sample2status.tsv raw.ped
do
	if test ! -s "\${F}"
	then
		rm -fv "\${F}"
	fi
done

touch versions.yml
"""
}
