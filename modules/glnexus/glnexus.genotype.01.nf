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

include {moduleLoad} from '../utils/functions.nf'


process GLNEXUS_GENOTYPE_01  {
tag "${bed.name} N=${vcfs.size()}"
memory "15g"
cpus 5
afterScript  "rm -rf TMP GLnexus.DB"
input:
	val(meta)
	tuple path(bed),val(vcfs)
output:
	tuple path(bed),path("merged.bcf"),emit:output
	path("version.xml"),emit:version
	path("merged.bcf.csi")
script:
       	def img = meta.glnexus_singularity_image?:"/LAB-DATA/BiRD/users/lindenbaum-p/packages/glnexus/glnexus.simg"
	def config = meta.glnexus_config?:"DeepVariantWGS"
"""
	hostname 1>&2
	${moduleLoad("bcftools")}
	mkdir TMP

cat << EOF > TMP/jeter.txt
${vcfs.join("\n")}
EOF

i=1
cat TMP/jeter.txt | while read F
do
	ln -vs "\${F}" "TMP/tmp.\${i}.vcf.gz"
	ln -vs "\${F}.tbi" "TMP/tmp.\${i}.vcf.gz.tbi"
	echo "/tmpdir/tmp.\${i}.vcf.gz" >> TMP/vcf.list
	((i=i+1))
	
done

       	singularity exec\
               	--home \${PWD} \
		--bind \${PWD}/TMP:/tmpdir \
               	--bind \${PWD}:/beddir \
               	--bind \${PWD}:/outdir \
               	${img} \
		/usr/local/bin/glnexus_cli \
		--config ${config} \
               	--bed /beddir/${bed.name} \
		${task.cpus?" --threads ${task.cpus}":""} \
		--mem-gbytes ${task.memory.giga} \
		--list /tmpdir/vcf.list > TMP/jeter.vcf.gz
	
	bcftools +fill-tags ${task.cpus?" --threads ${task.cpus}":""} -O b -o merged.bcf TMP/jeter.vcf.gz -- -t AN,AC,AF
	bcftools index ${task.cpus?" --threads ${task.cpus}":""} "merged.bcf"


##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">Combine GVCF with glnexus</entry>
        <entry key="bed">${bed}</entry>
        <entry key="gvcfs.count">${vcfs.size()}</entry>
   	<entry key="config">${config}</entry>
        <entry key="glnexus.version">\$( bcftools view --header-only merged.bcf | grep -m1 -o GLnexusVersion=.*) </entry>
</properties>
EOF

"""
}

