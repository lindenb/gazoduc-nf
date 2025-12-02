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
include { verify  } from '../../../modules/utils/functions.nf'
include { isBlank } from '../../../modules/utils/functions.nf'

workflow META_TO_PED{
take:
	metadata
	metas
main:
	versions = Channel.empty()
	multiqc = Channel.empty()
	MAKE_PED(
		metadata,
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
	versions = versions.mix(MAKE_PED.out.versions)
	
	cases = MAKE_PED.out.cases.ifEmpty([[id:"no.cases"],[]])
	controls = MAKE_PED.out.controls.ifEmpty([[id:"no.controls"],[]])
	pedigree = MAKE_PED.out.ped.ifEmpty([[id:"nojvarkitped"],[]])
	sample2collection = MAKE_PED.out.sample2col.ifEmpty( [[id:"nosn2col"],[]])
	pedigree_gatk  = MAKE_PED.out.gatk.ifEmpty( [[id:"nogatkped"],[]])
	sample2status = MAKE_PED.out.sample2status.ifEmpty( [[id:"nosample2status"],[]])
	males = MAKE_PED.out.males.ifEmpty([[id:"no.males"],[]])
	females = MAKE_PED.out.females.ifEmpty([[id:"no.females"],[]])
	
	MAKE_DEFAULT_SELECT_CODE(cases,controls)
	versions = versions.mix(MAKE_DEFAULT_SELECT_CODE.out.versions)
	select_code = MAKE_DEFAULT_SELECT_CODE.out.select_code.ifEmpty([[id:"no.code"],[]])
emit:
	cases
	controls
	pedigree
	pedigree_gatk
	sample2collection
	sample2status
	males
	females
	select_code
	versions
	multiqc
}

process MAKE_PED {
label "process_single"
conda "${moduleDir}/../../../conda/multiqc.yml"
input:
	val(meta)
	val(L)
output:
	tuple val(meta),path("*.pedigree4gatk.ped"),optional:true,emit:gatk
	tuple val(meta),path("*.jvarkit.ped"),optional:true,emit:ped
	tuple val(meta),path("*.males.txt"),optional:true,emit:males
	tuple val(meta),path("*.females.txt"),optional:true,emit:females
	tuple val(meta),path("*.cases.txt"),optional:true,emit:cases
	tuple val(meta),path("*.controls.txt"),optional:true,emit:controls
	tuple val(meta),path("*.sample2collection.tsv"),optional:true,emit:sample2col
	tuple val(meta),path("*.sample2status.tsv"),optional:true,emit:sample2status
	
	path("versions.yml"),emit:versions
script:

	def prefix= task.ext.prefix?:"${meta.id?:""}"
	verify(prefix!=null && !isBlank(prefix.toString()),"${task.process} meta.id must be defined")
"""
cat << EOF | sort -T . | uniq > raw.ped
${L.join("\n")}
EOF

python3 ${moduleDir}/ped.py raw.ped > /dev/null

for F in males.txt females.txt cases.txt controls.txt sample2collection.tsv pedigree4gatk.ped sample2status.tsv jvarkit.ped
do
	if test ! -s "\${F}"
	then
		rm -fv "\${F}"
	else
		mv -v "\${F}" "${prefix}.\${F}"
	fi
done

cat << EOF > versions.yml
${task.process}:
	python: todo
EOF
"""

stub:
	def prefix= task.ext.prefix?:"${meta.id?:""}"
	verify(prefix!=null && !isBlank(prefix.toString()),"${task.process} meta.id must be defined")
"""
for F in males.txt females.txt cases.txt controls.txt sample2collection.tsv pedigree4gatk.ped sample2status.tsv jvarkit.ped
do
	touch "${prefix}.\${F}"
done
touch versions.yml
"""
}


process MAKE_DEFAULT_SELECT_CODE {
executor "local"
input:
	tuple val(meta1),path(cases)
	tuple val(meta2),path(ctrls)
output:
	tuple val(meta1),path("default.select.code"),emit:select_code
	path("versions.yml"),emit:versions
script:
"""
echo 'return true;' > default.select.code

if test -s "${cases?:"NO_CASES"}" && test -s ${ctrls?:"NO_CTRL"}
then

	echo 'final Set<String> cases = new HashSet<>(Arrays.asList(' > jeter.code

	sed 's/^/"/;s/\$/"/' '${cases}' | paste -sd ',' >> jeter.code

	echo '));' >> jeter.code


	echo 'final Set<String> controls = new HashSet<>(Arrays.asList(' >> jeter.code

	sed 's/^/"/;s/\$/"/' '${ctrls}' | paste -sd ',' >> jeter.code

	echo '));' >> jeter.code

cat << '__EOF__' >>  jeter.code
for(Genotype g: variant.getGenotypes()) {
        final boolean has_alt = g.hasAltAllele();
        if(!has_alt &&  cases.contains(g.getSampleName())) return false;
        if( has_alt &&  controls.contains(g.getSampleName())) return false;
        }
return true;
__EOF__

mv jeter.code default.select.code

fi


touch versions.yml
"""
stub:
"""
echo 'return true;' > default.select.code
touch versions.yml
"""
}
