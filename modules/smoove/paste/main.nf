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
process SMOOVE_PASTE {
tag "${meta.id}"
label "process_short"
label "smoove"
afterScript  "rm -rf TMP"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta ),path("VCFS/*")
output:
	tuple val(meta ),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:meta.id+".paste"
	def limit = (task.ext.limit?:200) as int
"""
	hostname 1>&2
	# smoove will write to the system TMPDIR. For large cohorts, make sure to set this to something with a lot of space
	mkdir -p TMP/TMP2
	export TMPDIR=\${PWD}/TMP/TMP2
	set -x

	find VCFS/ -name "*.vcf.gz" | sort -T TMP > TMP/jeter.list

	if test  `wc -l < TMP/jeter.list` -le ${limit}
	then

		smoove paste --name "${prefix}" VCFS/*.vcf.gz

	else
	
		SQRT=`awk 'END{X=NR;if(X<10){print(X);} else {z=sqrt(X); print (z==int(z)?z:int(z)+1);}}' TMP/jeter.list`
		split -a 9 --additional-suffix=.list --lines=\${SQRT} TMP/jeter.list TMP/chunck.
		i=1

		find TMP/ -type f -name "chunck*.list" | while read F
		do

			smoove paste --name "jeter" \$(cat "\${F}")
			bcftools index --threads ${task.cpus} -t -f "jeter.smoove.square.vcf.gz"

			mv "jeter.smoove.square.vcf.gz" TMP/\${i}.square.vcf.gz
			mv "jeter.smoove.square.vcf.gz.tbi" TMP/\${i}.square.vcf.gz.tbi
			
			echo "TMP/\${i}.square.vcf.gz" >> TMP/jeter2.list
			i=\$((i+1))		
		done
		
		smoove paste --name "${prefix}" \$(cat  TMP/jeter2.list)

	fi

bcftools sort  -T TMP/sort  --max-mem "${task.memory.giga}G" -O z -o TMP/jeter.vcf.gz  "${prefix}.smoove.square.vcf.gz"
mv -v TMP/jeter.vcf.gz "${prefix}.smoove.square.vcf.gz"
bcftools index --threads ${task.cpus} -t -f "${prefix}.smoove.square.vcf.gz"


cat <<- EOF > versions.yml
"${task.process}":
    smoove: \$(smoove 2>&1 | awk '(NR==1) {print \$3;}')
EOF
"""

stub:
def prefix = task.ext.prefix?:meta.id+".paste"
"""
touch versions.yml ${prefix}.smoove.square.vcf.gz ${prefix}.smoove.square.vcf.gz.tbi
"""
}
