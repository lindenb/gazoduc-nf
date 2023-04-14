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


include {moduleLoad;getVersionCmd} from './../../modules/utils/functions.nf'
include {gatkGetArgumentsForCombineGVCFs;gatkGetArgumentsForGenotypeGVCF} from './gatk.hc.utils.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'

workflow GATK_PARALLEL_COMBINE_GVCFS {
	take:
		meta
		reference
		pedigree
		rows
	main:
		version_ch = Channel.empty()


		level0_ch = COMBINE_LEVEL0(meta, reference, ch1.output )
		version_ch = version_ch.mix(level0_ch.version)

		level1_ch = COMBINE_LEVEL1(meta, reference, level0_ch.output.groupTuple())
		version_ch = version_ch.mix(level1_ch.version)

		level2_ch = GENOTYPE_LEVEL2(meta, reference, pedigree, level1_ch.output )
		version_ch = version_ch.mix(level1_ch.version)

		version_ch = MERGE_VERSION(meta, "combine_gvcfs", "combine.gvcfs", version_ch.collect())

	emit:
		region_vcf = level2_ch.output
		version = version_ch
	}


process COMBINE_LEVEL0 {
tag "${row.interval} ${row.gvcf_split}"
memory {task.attempt <2 ? "15g":"60g"}
memory "10g"
cpus "3"
errorStrategy  'retry'
maxRetries 3
afterScript 'rm -rf  TMP'
input:
	val(meta)
	val(reference)
	val(row)
output:
        tuple val("${row.interval}"),path("combine0.g.vcf.gz"),emit:output
        path("combine0.g.vcf.gz"),emit:index
	path("version.xml"),emit:version
script:
     def otherOpts =  gatkGetArgumentsForCombineGVCFs(meta.plus("reference":reference))
"""
hostname 1>&2
${moduleLoad("gatk4")}

mkdir -p TMP

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" CombineGVCFs \
                ${otherOpts} \
		-L "${row.interval}" \
		-V "${row.gvcf_split}"  \
		-O "TMP/combine0.g.vcf.gz"

mv -v TMP/combine0.g.vcf.gz ./
mv -v TMP/combine0.g.vcf.gz.tbi ./

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">combine level0</entry>
        <entry key="region">${row.interval}</entry>
	<entry key="versions">${getVersionCmd("gatk")}</entry>
</properties>
EOF
"""
}


process COMBINE_LEVEL1 {
tag "${row.interval} N=${L.size()}"
memory {task.attempt <2 ? "15g":"60g"}
memory "10g"
cpus "3"
errorStrategy  'retry'
maxRetries 3
afterScript 'rm -rf  TMP'
input:
	val(meta)
	val(reference)
	tuple val(region),val(L)
output:
        tuple val("${row.interval}"),path("combine0.g.vcf.gz"),emit:output
        path("combine0.g.vcf.gz"),emit:index
	path("version.xml"),emit:version
script:
     def otherOpts =  gatkGetArgumentsForCombineGVCFs(meta.plus("reference":reference)) 
"""
hostname 1>&2
${moduleLoad("gatk4")}
mkdir -p TMP

if [ "\${L.size()}" == "1" ] ; then

	ln -s "${L[0].toRealPath()}" combine1.g.vcf.gz
	ln -s "${L[0].toRealPath()}.tbi" combine1.g.vcf.gz.tbi

	sleep 5
else

cat << EOF > TMP/vcfs.list
${L.join("\n")}
EOF

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" CombineGVCFs \
        ${otherOpts} \
        -L "${region}" \
	-V TMP/vcfs.list  \
	-O "TMP/combine1.g.vcf.gz"

mv TMP/combine1.g.vcf.gz ./
mv TMP/combine1.g.vcf.gz.tbi ./

fi

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">combine level 1</entry>
        <entry key="region">${row.interval}</entry>
	<entry key="versions">${getVersionCmd("gatk")}</entry>
</properties>
EOF
"""
}

process GENOTYPE_LEVEL2 {
tag "${region}"
memory {task.attempt <2 ? "15g":"60g"}
memory "10g"
cpus "3"
errorStrategy  'retry'
maxRetries 3
afterScript 'rm -rf  TMP'
input:
        val(meta)
        val(reference)
        path(pedigree)
        tuple val(region),path(gvcf)
output:
        tuple val("${region}"),path("genotyped.bcf"),emit:output
        path("genotyped.bcf.csi"),emit:index
        path("version.xml"),emit:version
script:
     def otherOpts =  gatkGetArgumentsForGenotypeGVCF(meta.plus("reference":reference,"pedigree":pedigree))

"""
hostname 1>&2
${moduleLoad("gatk4 bcftools")}
mkdir -p TMP

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" GenotypeGVCFs  \
      ${otherOpts} \
      -L "${region}" \
      -V "${gvcf.toRealPath()}" \
      -O "TMP/jeter.vcf.gz"

bcftools view --threads ${task.cpus} -O b -o genotyped.bcf TMP/jeter.vcf.gz
bcftools index --threads ${task.cpus} genotyped.bcf

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">genotype gvcfs</entry>
        <entry key="region">${region}</entry>
	<entry key="versions">${getVersionCmd("gatk bcftools")}</entry>
</properties>
EOF
"""
}