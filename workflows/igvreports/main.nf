/*

Copyright (c) 2026 Pierre Lindenbaum

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
nextflow.enable.dsl=2


include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {dumpParams;runOnComplete;moduleLoad} from '../../modules/utils/functions.nf'
//include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {DOWNLOAD_CYTOBAND} from '../../modules/ucsc/download.cytoband.nf'
include {DOWNLOAD_REFGENE } from '../../modules/ucsc/download.refgene.nf'

if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}



workflow {
	
	if(params.fasta==null) {
      throw new IllegalArgumentException("undefined --fasta");
      }
	  
	def genome_hash= [
		id:file(params.fasta).simpleName,
		name:file(params.fasta).simpleName
		]

	ch1 = IGVREPORTS(
		[id:"igvreports"],
		[genome_hash,file(params.fasta)],
		[genome_hash,file(params.fai)],
		[genome_hash,file(params.dict)],
		file(params.bams),
		file(params.vcf)
		)
	//html = VERSION_TO_HTML(ch1.version)
	}

runOnComplete(workflow)

workflow IGVREPORTS {
	take:
		meta
		fasta
		fai
		dict
		bams
		vcf
	main:
		version_ch = Channel.empty()


        snbam_ch = SAMTOOLS_SAMPLES(fasta,fai,dict)
		

		compile_ch = COMPILE()
		version_ch = version_ch.mix(compile_ch.version)

		cyto_ch = DOWNLOAD_CYTOBAND(fasta,fai,dict)

		refgene_ch = DOWNLOAD_REFGENE(fasta,fai,dict)

		prepare_ch = PREPARE_IGVREPORTS(
			fasta,
			fai,
			dict
			compile_ch.output,
			snbam_ch.output, 
			vcf
			)
                version_ch = version_ch.mix(prepare_ch.version)

		report_ch = IGVREPORT(fasta,fai,dict, vcf, cyto_ch.output, refgene_ch.output, prepare_ch.output.splitCsv(header:true,sep:'\t') )
                version_ch = version_ch.mix(report_ch.version)
		
		ZIPIT(report_ch.output.collect())

	}

process SAMTOOLS_SAMPLES {
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
output:
	path("sample_bam.tsv"),emit:output
script:
"""
set -o pipefail
samtools samples -f '${fasta}' < "${bams}" | awk -F '\t' '\$3!="."' | cut -f1,2 > sample_bam.tsv
test -s sample_bam.tsv
"""
}

process COMPILE {
        afterScript "rm -rf TMP"
        output:
                path("minikit.jar"),emit:output
                path("version.xml"),emit:version
        script:
        """
        hostname 1>&2
        ${moduleLoad("jvarkit")}

cat << EOF > Minikit.java

EOF

mkdir -p TMP
javac -d TMP -cp \${JVARKIT_DIST}/vcffilterjdk.jar:. Minikit.java
jar cvf minikit.jar -C TMP .


###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">compile minikit</entry>
        <entry key="javac.version">\$(javac -version 2>&1)</entry>
</properties>
EOF
"""
}

process PREPARE_IGVREPORTS {
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	path(jar)
	path(sn2bam)
	path(vcf)
output:
	path("output.tsv"),emit:output
	path("version.xml"),emit:version
script:
"""        
hostname 1>&2     
${moduleLoad("jvarkit bcftools")}
set -o pipefail
mkdir -p OUT



cat ${pricess } sed 's/__NUM_CONTROLS__/${params.num_controls/

bcftools view "${vcf}" |\\
	java -cp ${jar}:\${JVARKIT_DIST}/jvarkit.jar Minikit \\
			--reference "${fasta}" \\
			--bams ${sn2bam} \\
			--out "\${PWD}/OUT" > output.tsv


###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">prepare data for IGV reports</entry>
        <entry key="versions"></entry>
</properties>
EOF
"""
}

process IGVREPORT {
tag "${row.title}"
afterScript "rm -rf TMP"
conda "${params.conda}/IGVREPORTS"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	path(vcf)
	tuple val(meta5),path(cytoband)
	tuple val(meta6),path(refgene)
	val(row)
output:
	path("${row.title}.html"),emit:output
	path("version.xml"),emit:version
script:
	def refgene2 = refgene.find{it.name.endsWith(".txt.gz")}.join(" ")
"""
hostname 1>&2
mkdir -p TMP

create_report ${row.vcf.isEmpty()?row.bedpe:row.vcf}  ${fasta} \\
	--ideogram "${cytoband}" \\
	${row.flanking.isEmpty()?"":"--flanking ${row.flanking}"} \\
	${row.info.isEmpty()?"":"--info-columns ${row.info}"} \\
	--tracks ${vcf} ${row.bams} ${refgene2} \\
	--output TMP/${row.title}.html

mv -v "TMP/${row.title}.html" ./

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">generate IGV reports</entry>
        <entry key="versions"></entry>
</properties>
EOF
"""
}

process ZIPIT {
input:
	path(htmls)
output:
	path("archive.zip"),emit:output
script:
"""
zip -9 -j archive.zip ${htmls}
"""
}
