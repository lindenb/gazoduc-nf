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

include { SAMTOOLS_SAMPLES01 as CASES_BAMS; SAMTOOLS_SAMPLES01 as CTRLS_BAMS} from '../../modules/samtools/samtools.samples.01.nf'

workflow SAMTOOLS_CASES_CONTROLS_01 {
	take:
		meta
		reference
		cases_bams
		ctrls_bams
	main:
		version_ch = Channel.empty()
		
		cases_ch = CASES_BAMS([:], reference, cases_bams)
		version_ch = version_ch.mix(cases_ch.version)

		if(ctrls_bams.isEmpty()) {
			merge_ch = MERGE_SN([:],cases_ch.out,"")
			}
		else
			{
			ctrls_ch = CTRLS_BAMS([:], reference, ctrls_bams)
			version_ch = version_ch.mix(ctrls_ch.version)
			merge_ch = MERGE_SN([:],cases_ch.out,ctrls_ch.out)
			}
	emit:
		output = merge_ch.out
		version = version_ch
	}

process MERGE_SN {
executor "local"
input:
	val(meta)
	val(cases_list)
	val(ctrls_list)
output:
	path("cases_controls_bams.tsv"),emit:out
	path("version.xml"),emit:version
script:
"""

awk '{printf("%s\tcase\\n",\$0);}' "${cases_list}" > jeter.tsv

if [ -f "${ctrls_list}" ] ; then

	awk '{printf("%s\tcontrol\\n",\$0);}' "${ctrls_list}" >> jeter.tsv

fi

# no dup samples
cut -f 1 jeter.tsv | sort -T . | uniq -d > jeter.txt
test ! -s jeter.txt

# no dup bam
cut -f 2 jeter.tsv | sort -T . | uniq -d > jeter.txt
test ! -s jeter.txt

sort -T . -t '\t' -k1,1 jeter.tsv > cases_controls_bams.tsv


rm jeter.tsv jeter.txt


#############################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">merge sample/bam list for cases and controls</entry>
</properties>
EOF
"""
}
