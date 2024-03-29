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


include {moduleLoad;getVersionCmd} from './../../modules/utils/functions.nf'
include {gatkGetArgumentsForCombineGVCFs;gatkGetArgumentsForGenotypeGVCF} from './gatk.hc.utils.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'

workflow GATK_PARALLEL_COMBINE_GVCFS {
	take:
		meta
		genomeId
		pedigree
		rows
	main:
		version_ch = Channel.empty()

		ch1 = EXPAND(meta,rows)		
		version_ch = version_ch.mix(ch1.version)

		level0_ch = COMBINE_LEVEL0(meta, genomeId, ch1.output.splitCsv(header:true,sep:'\t') )
		version_ch = version_ch.mix(level0_ch.version)

		level1_ch = COMBINE_LEVEL1(meta, genomeId, level0_ch.output.groupTuple())
		version_ch = version_ch.mix(level1_ch.version)

		level2_ch = GENOTYPE_LEVEL2(meta, genomeId, pedigree, level1_ch.output )
		version_ch = version_ch.mix(level1_ch.version)

		version_ch = MERGE_VERSION(meta, "combine_gvcfs", "combine.gvcfs", version_ch.collect())

	emit:
		region_vcf = level2_ch.output
		version = version_ch
	}


process EXPAND {
executor "local"
tag "${row.interval} ${row.gvcf_split}"
input:
	val(meta)
	val(row)
output:
	path("output.tsv"),emit:output
	path("version.xml"),emit:version
script:
"""
awk 'BEGIN{printf("interval\tgvcfs\\n");} {printf("${row.interval}\t%s\\n",\$0);}' "${row.gvcf_split}" > output.tsv

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">extract gvcf list</entry>
        <entry key="region">${row.interval}</entry>
	<entry key="versions">${getVersionCmd("awk")}</entry>
</properties>
EOF
"""
}

process COMBINE_LEVEL0 {
tag "${row.interval} ${row.gvcfs}"
memory {task.attempt <2 ? "15g":"60g"}
memory "10g"
cpus "1"
errorStrategy  'retry'
maxRetries 3
afterScript 'rm -rf  TMP'
input:
	val(meta)
	val(genomeId)
	val(row)
output:
        tuple val("${row.interval}"),path("combine0.g.vcf.gz"),emit:output
        path("combine0.g.vcf.gz"),emit:index
	path("version.xml"),emit:version
script:
     def otherOpts =  gatkGetArgumentsForCombineGVCFs(meta.plus("genomeId":genomeId))
"""
hostname 1>&2
${moduleLoad("gatk4")}

mkdir -p TMP

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" CombineGVCFs \
                ${otherOpts} \
		-L "${row.interval}" \
		-V "${row.gvcfs}"  \
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
tag "${region} N=${L.size()}"
memory {task.attempt <2 ? "15g":"60g"}
memory "10g"
cpus "3"
errorStrategy  'retry'
maxRetries 3
afterScript 'rm -rf  TMP'
input:
	val(meta)
	val(genomeId)
	tuple val(region),val(L)
output:
        tuple val("${region}"),path("combine1.g.vcf.gz"),emit:output
        path("combine1.g.vcf.gz"),emit:index
	path("version.xml"),emit:version
script:
     def otherOpts =  gatkGetArgumentsForCombineGVCFs(meta.plus("genomeId":genomeId)) 
if(L.size()==1)
"""

ln -s "${L[0].toRealPath()}" combine1.g.vcf.gz
ln -s "${L[0].toRealPath()}.tbi" combine1.g.vcf.gz.tbi

echo "<properties/>" > version.xml
"""
else
"""
hostname 1>&2
${moduleLoad("gatk4")}
mkdir -p TMP

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


##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">combine level 1</entry>
        <entry key="region">${region}</entry>
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
        val(genomeId)
        path(pedigree)
        tuple val(region),path(gvcf)
output:
        tuple val("${region}"),path("genotyped.bcf"),emit:output
        path("genotyped.bcf.csi"),emit:index
        path("version.xml"),emit:version
script:
     def otherOpts =  gatkGetArgumentsForGenotypeGVCF(meta.plus("genomeId":genomeId,"pedigree":pedigree))

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
