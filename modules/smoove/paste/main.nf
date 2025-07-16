
process SMOOVE_PASTE {
tag "${meta.id}"
label "process_short"
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
	def limit = 200
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


bcftools index --threads ${task.cpus} -t -f "${prefix}.smoove.square.vcf.gz"


cat <<- EOF > versions.yml
"${task.process}":
    smoove: todo
EOF
"""
}
