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
include {DUPHOLD as RUN_DUPHOLD  } from '../../modules/duphold'
include {BCFTOOLS_MERGE          } from '../../modules/bcftools/merge'
/*
include {SAMTOOLS_SAMPLES_01} from '../samtools/samtools.samples.01.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {moduleLoad;isHg19;isHg38;runOnComplete;parseBoolean;getVersionCmd} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
*/

workflow DUPHOLD {
	take:
		meta
		fasta
		fai
		dict
		bams //meta,bam,bai,[fasta,fai,dict]?
		vcf /* SV/CNV vcf */
		snps /** small indels , snps */
	main:
		versions  = Channel.empty()
		multiqc = Channel.empty()

		ch1 =  bams.branch{v->
			without_fasta: v.size()==3 //meta,bam,bai
			with_fasta: v.size()==6 //meta,bam,bai,fasta,fai,dict
			others: true
			}
		
		ch1.others.view{"BOUM CHECK DUPHOLD: ${it}."}

		ch2 = ch1.without_fasta.combine(fasta).combine(fai) // meta,bam,bai,meta,fasta,meta,fai,meta,dict
			.map{[it[0],it[1],it[2],it[4],it[6]]}
			.mix(ch1.with_fasta.map{[it[0],it[1],it[2],it[3],it[4]]})
			.combine(vcf.map{[it[1]]})
			.combine(snps.map{[it[1],it[2]]})
			.multiMap{
				fasta: [[id:it[3].name],it[3]]
				fai: [[id:it[3].name],it[4]]
				args: [it[0],it[1],it[2],it[5],it[6],it[7]]
			}
	ch2.args.view()
		RUN_DUPHOLD(
			ch2.fasta,
			ch2.fai,
			ch2.args
			)
	versions = versions.mix(RUN_DUPHOLD.out.versions)

	BCFTOOLS_MERGE(RUN_DUPHOLD.out.vcf
		.map{[it[1],it[2]]}
		.collec()
		.map{[[id:"duphold"],it]}
		)
	versions = versions.mix(BCFTOOLS_MERGE.out.versions)
	
	emit:
		versions
		vcf = BCFTOOLS_MERGE.out.vcf
	}


process DISPATCH_VCF {
tag "${contig} ${file(vcf).name}"
input:
	val(meta)
	tuple val(contig),val(vcf)
output:
	path("vcfs.list"),emit:output
	path("remain.${contig}.bcf"),emit:remain
	path("version.xml"),emit:version
script:
	def method = meta.duphold_split_method?:"--variants-count 100"
"""
hostname 1>&2
${moduleLoad("bcftools jvarkit")}
set -o pipefail
mkdir -p OUT

bcftools view -i 'SVTYPE=="DEL" || SVTYPE=="DUP"' "${vcf}" "${contig}" |\
	java -jar \${JVARKIT_DIST}/vcfsplitnvariants.jar ${method} -o \${PWD}/OUT/split.${contig}.

find OUT/ -type f -name "*.vcf.gz" | while read F
do
	bcftools index -t --force "\${F}"
done

find \${PWD}/OUT/ -type f -name "*.vcf.gz" > vcfs.list

# remaining
bcftools view -e 'SVTYPE=="DEL" || SVTYPE=="DUP"' "${vcf}" "${contig}" -O b -o "remain.${contig}.bcf"
bcftools index "remain.${contig}.bcf"

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">split VCF per contig and using jvarkit/vcfsplitnvariants</entry>
        <entry key="method">${method}</entry>
        <entry key="vcf">${vcf}</entry>
        <entry key="contig">${contig}</entry>
        <entry key="versions">${getVersionCmd("bcftools jvarkit/vcfsplitnvariants")}</entry>
</properties>
EOF
"""
}

process JOIN_BAM_VCFS {
executor "local"
tag "${samples} ${bams}"
afterScript "rm -rf TMP"
input:
	val(meta)
	path(bams)
	path(samples)
output:
	path("join.bams.list"),emit:bams
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
mkdir -p TMP

cut -f1,3 "${bams}" | sort -t '\t' -k1,1 -T TMP > TMP/jeter.a

sort -t '\t' -k1,1 -T TMP "${samples}"> TMP/jeter.b

join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter.a TMP/jeter.b > join.bams.list

test -s join.bams.list

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">join vcf samples and bams</entry>
        <entry key="samples">${samples}</entry>
        <entry key="bams">${bams}</entry>
</properties>
EOF
"""
}

