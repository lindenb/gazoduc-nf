/*

Copyright (c) 2024 Pierre Lindenbaum

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
include {moduleLoad;getVersionCmd;parseBoolean} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {JVARKIT_VCF_TO_INTERVALS_01} from './jvarkit.vcf2intervals.nf'




/**
 * call jvarkit vcfgatkeval for one or more vcf
 *
 */
workflow JVARKIT_GATK_HARD_FILTERING_02 {
	take:
		reference
		vcf /* path to one or more vcfs */
		bed /* limit to that BED or NO_FILE */
	main:
		if(!meta.containsKey("percentile")) throw new IllegalArgumentException("percentile undefined");
		if(!meta.containsKey("hard_filtering")) throw new IllegalArgumentException("hard_filtering");
		
		version_ch = Channel.empty()
		intervals_ch = JVARKIT_VCF_TO_INTERVALS_01(vcf, bed)
		version_ch = version_ch.mix(intervals_ch.version)

		contig_vcf_ch = intervals_ch.bed.splitCsv(sep:'\t',header:false)

		rch = FOR_EACH_INTERVAL(meta, contig_vcf_ch)
		version_ch = version_ch.mix(rch.version)

		concat_ch = CONCAT_TABLES(rch.output.collect())
		version_ch = version_ch.mix(concat_ch.version)

		flt_ch = APPLY_FILTERS(genomeId, contig_vcf_ch.combine(concat_ch.output))
		version_ch = version_ch.mix(flt_ch.version)

		concat_ch = CONCAT(flt_ch.bed.collect())


		version_ch = MERGE_VERSION("gatkeval",version_ch.collect())
	emit:
		bed = concat_ch.bed
		pdf = concat_ch.pdf
		version = version_ch
	}

process FOR_EACH_INTERVAL {
tag "${contig}:${start}-${end} ${file(vcf).name}"
afterScript "rm -rf TMP"
memory '2g'
input:
	val(meta)
	tuple val(contig),val(start),val(end),val(vcf)
output:
	path("OUT/gatk.eval.output.table.txt"),emit:output
	path("version.xml"),emit:version
script:
	if(!meta.containsKey("percentile")) throw new IllegalArgumentException("percentile undefined");
	def percentil = meta.percentil
"""
hostname 1>&2
${moduleLoad("jvarkit bcftools")}
set -o pipefail
mkdir -p TMP

echo '${chrom}\t${start}\t${end}' > TMP/jeter.bed

bcftools view -G --regions-file TMP/jeter.bed '${vcf}' |\
	java -jar -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP \${JVARKIT_DIST}/jvarkit.jar vcfgatkeval --percentile ${meta.percentile} --input-type vcf -o TMP/gatk.eval
	
mv TMP OUT
	
###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="versions">${getVersionCmd("bcftools jvarkit/jvarkit")}</entry>
</properties>
EOF
"""
}

process CONCAT_TABLES {
tag "N=${L.size()}"
afterScript "rm -rf TMP"
memory '2g'
input:
	val(L)
output:
	path("OUT/gatk.eval.output.filters.txt"),emit:output
	path("OUT/gatk.eval.output.pdf"),emit:pdf
	path("version.xml"),emit:version	
script:
"""
hostname 1>&2
${moduleLoad("jvarkit R")}
set -o pipefail
mkdir -p TMP

cat << EOF > TMP/jeter.list
${L.join("\n")}
EOF

xargs -a TMP/jeter.list -L 50 cat | java -jar -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP \${JVARKIT_DIST}/jvarkit.jar vcfgatkeval --percentile ${meta.percentile} --input-type table -o TMP/gatk.eval

sed 's%output.pdf%TMP/gatk.eval.output.pdf%' TMP/gatk.eval.output.R | R --vanilla 

mv TMP OUT
###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="versions">${getVersionCmd("bcftools jvarkit/jvarkit")}</entry>
</properties>
EOF
"""
}


process APPLY_FILTERS {
tag "${contig}:${start}-${end} ${file(vcf).name}"
afterScript "rm -rf TMP"
memory '3g'
input:
	val(genomeId)
	tuple val(contig),val(start),val(end),val(vcf),path(filters)
output:
	path("variants.bed"),emit:bed
	path("version.xml"),emit:version	
script:
	if(!meta.containsKey("hard_filtering")) throw new IllegalArgumentException("hard_filtering undefined");
	def reference = params.genomes[genomeId]
"""
hostname 1>&2
${moduleLoad("gatk/0.0.0 bcftools")}
set -o pipefail
mkdir -p TMP

echo '${chrom}\t${start}\t${end}' > TMP/jeter.bed

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" VariantFiltration \
	-L TMP/jeter.bed \
	-V '${vcf}' \
	-R '${reference}' \
	-O TMP/variants.vcf.gz \
	--arguments_file ${filters}

if ${parseBoolean(meta.hard_filtering)} ; then
	bcftools view --apply-filters 'PASS,.' -O z -o TMP/variants2.vcf.gz TMP/variants.vcf.gz
	mv -v TMP/variants2.vcf.gz TMP/variants.vcf.gz
	bcftools index -t -f TMP/variants.vcf.gz
fi

mv TMP/variants.vcf.gz ./
mv TMP/variants.vcf.gz.tbi ./
touch variants.vcf.gz.tbi

echo "${chrom}\t${start}\t${end}\t\${PWD}/variants.vcf.gz" > variants.bed

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
</properties>
EOF
"""
}

process CONCAT {
executor "local"
afterScript "rm -rf TMP"
tag "N=${L.size()}"
input:
	val(L)
output:
	path("variants.bed"),emit:bed
	path("vcfs.list"),emit:vcfs
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP

cat << EOF > TMP/jeter.list
${L.join("\n")}
EOF

xargs -a TMP/jeter.list -L 50 cat | sort -T TMP -t '\t' -k1,1 -k2,2n | uniq > TMP/variants.bed
mv TMP/variants.bed ./

cut -f 4 variants.bed | sort | uniq > vcfs.list

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
</properties>
EOF
"""	
}
