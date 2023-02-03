/*

Copyright (c) 2023 Pierre Lindenbaum

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

def gazoduc = gazoduc.Gazoduc.getInstance()

gazoduc.make("manta_cpu",16).
        description("Number of cpus for manta").
        setInt().
        put()


include { getKeyValue; getModules} from '../../modules/utils/functions.nf'

process MANTA_SINGLE_01 {
    tag "${name} ${file(bam).name}"
    afterScript "rm -rf TMP"
    cache 'lenient'
    errorStrategy 'finish'
    cpus ((params.manta_cpus?:16) as int)
    memory "20g"
    input:
	val(meta)
	val(reference)
	tuple val(name),val(bam)
    output:
    	tuple val(name),val(bam),path("${name}.txt"),emit:output
	path("version.xml"),emit:version
    script:
	def prefix = meta.getOrDefault("prefix","")
	"""
	hostname 1>&2
	module load ${getModules("manta samtools bcftools")}
	mkdir -p TMP

	configManta.py  --bam "${bam}" --referenceFasta "${reference}" \
		--runDir "TMP"

	
	./TMP/runWorkflow.py --quiet -m local -j ${task.cpus}
	
	rm -rf ./TMP/workspace


	# convert BND TO INVERSIONS (added 20230115 but not tested)
	DIPLOID=`find ./TMP -type f -name "*.diploidSV.vcf.gz"`
	test ! -z "\${DIPLOID}"
	\$(ls \$( dirname \$(which configManta.py) )/../share/manta*/libexec/convertInversion.py)  `which samtools` "${reference}" "\${DIPLOID}" | bcftools sort -T TMP -O z -o TMP/jeter.vcf.gz

	bcftools index -t TMP/jeter.vcf.gz

	mv -v TMP/jeter.vcf.gz "\${DIPLOID}"
	mv -v TMP/jeter.vcf.gz.tbi "\${DIPLOID}.tbi"

	# change name to sample
	find ./TMP -type f -name "*.vcf.gz" \
		-printf 'mv -v %p  ${prefix}${name}.%f\\n mv -v %p.tbi ${prefix}${name}.%f.tbi\\n' |bash 
	
	ls *.vcf.gz
	find \${PWD}  -maxdepth 1  -type f -name "*.vcf.gz" -o -name "*.vcf.gz.tbi" |\
		awk '{printf("${name}\t${bam}\t%s\\n",\$0);}' > ${name}.txt
	grep -m1 vcf ${name}.txt



#################################################################################################

cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="date">\$(date)</entry>
	<entry key="sample">${name}</entry>
	<entry key="bam">${bam}</entry>
	<entry key="manta.version">\$(configManta.py --version)</entry>
</properties>
EOF
	"""
stub:
	"""
	touch "${name}.txt"
	echo "<properties/>" > version.xml
	"""
	}
