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
include {moduleLoad;getKeyValue;getModules;getVersionCmd} from '../../modules/utils/functions.nf'
include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'


workflow BCFTOOLS_CONCAT_PER_CONTIG_01 {
	take:
		/* params */
		meta
		/* row.vcfs a FILE containing the path to the indexed VCF */
		vcfs
	main:
		if(!meta.containsKey("suffix")) throw new IllegalArgumentException("meta.suffix is missing");
		if(!(meta.suffix.equals(".vcf.gz") || meta.suffix.equals(".bcf"))) throw new IllegalArgumentException("bad VCF suffix ${suffix}");
		
		version_ch = Channel.empty()

		c1_ch = VCF_PER_CONTIG(meta,vcfs)
		version_ch = version_ch.mix(c1_ch.version)

		
		c2_ch = c1_ch.output.splitCsv(header:false,sep:"\t").map{T->["contig":T[0],"vcfs":T[1]]}

		c3_ch = CONCAT_ONE_CONTIG(meta,c2_ch)
		version_ch = version_ch.mix(c3_ch.version)

		file_list_ch = COLLECT_TO_FILE_01([:],c3_ch.vcf.map{T->T.toRealPath()}.collect())

		c4_ch = VCF_TO_BED([with_header:false],file_list_ch.output)
		version_ch = version_ch.mix(c4_ch.version)

                version_ch = MERGE_VERSION("concat / contig", version_ch.collect())

	emit:
		version = version_ch
		vcfs = file_list_ch.output
		bed = c4_ch.bed		
	}




process VCF_PER_CONTIG {
tag "${vcfs}"
input:
	val(meta)
	path(vcfs)
output:
	path("ctg2vcfs.tsv"),emit:output	
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}
set -o pipefail

cat "${vcfs}" | while read V
do
	bcftools index -s "\${V}" | awk -F '\t' -vV=\${V} '{printf("%s\t%s\\n",\$1,V);}' >> ctg_vcf.txt
done

cut -f 1 ctg_vcf.txt | sort | uniq | while read C
do
	awk -vC=\$C -F '\t' '(\$1==C)' ctg_vcf.txt | cut -f2 |\
		sort -T . | uniq >> "contig.\${C}.vcf.list"

	echo "\$C\t\${PWD}/contig.\${C}.vcf.list" >> ctg2vcfs.tsv
done

rm ctg_vcf.txt
touch -c ctg2vcfs.tsv

###############################################

cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
        <entry key="description">build list of VCFs per contig</entry>
        <entry key="count(vcfs.in)">\$(wc -l < "${vcfs}")</entry>
        <entry key="count(contigs.out)">\$(wc -l < ctg2vcfs.tsv)</entry>
        <entry key="contigs">\$(cut -f 1 ctg2vcfs.tsv | paste -s -d,)</entry>
        <entry key="versions">${getVersionCmd("bcftools awk")}</entry>
</properties>
EOF

"""
}

process CONCAT_ONE_CONTIG {
tag "${row.contig}:${file(row.vcfs).name}"
cpus 3
input:
        val(meta)
        val(row)
output:
        path("${meta.prefix?:""}${row.contig}.concat${meta.suffix}"),emit:vcf
        path("${meta.prefix?:""}${row.contig}.concat${meta.suffix}${meta.suffix.contains("b")?".csi":".tbi"}"),emit:index
	path("version.xml"),emit:version
script:
	def prefix = meta.prefix?:""
	def suffix = meta.suffix?:".bcf"
	def contig = row.contig
	def vcfs = row.vcfs
	if(!(suffix.equals(".vcf.gz") || suffix.equals(".bcf"))) throw new IllegalArgumentException("bad VCF suffix ${suffix}");
        """
        hostname 1>&2
        ${moduleLoad("bcftools")}

	test -s "${vcfs}"

        bcftools concat --threads "${task.cpus}" --regions "${contig}" \
                --no-version --allow-overlaps --remove-duplicates \
                -O "${suffix.contains("b")?"b":"z"}" -o "${prefix}${contig}.concat${suffix}" --file-list "${vcfs}"

        bcftools index --threads "${task.cpus}" ${!suffix.contains("b")?"--tbi":""} "${prefix}${contig}.concat${suffix}"

###############################################

cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
        <entry key="description">concatenate VCF for one contig</entry>
        <entry key="count(vcfs.in)">\$(wc -l < "${vcfs}")</entry>
        <entry key="count(contigs)">${contig}</entry>
        <entry key="versions">${getVersionCmd("bcftools")}</entry>
</properties>
EOF

        """  
        }
