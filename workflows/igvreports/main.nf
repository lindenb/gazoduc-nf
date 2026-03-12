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

include { validateParameters                        } from 'plugin/nf-schema'
include { paramsHelp                                } from 'plugin/nf-schema'
include { paramsSummaryLog                          } from 'plugin/nf-schema'
include { parseBoolean                              } from '../../modules/utils/functions.nf'
include { assertKeyExistsAndNotEmpty                } from '../../modules/utils/functions.nf'
include { READ_SAMPLESHEET  as READ_BAMS            } from '../../subworkflows/nf/read_samplesheet'
include { READ_SAMPLESHEET  as READ_VARIANTS        } from '../../subworkflows/nf/read_samplesheet'
include { PREPARE_ONE_REFERENCE                     } from '../../subworkflows/samtools/prepare.one.ref'
include { META_TO_BAMS                              } from '../../subworkflows/samtools/meta2bams1'
include { IGV_REPORT                                } from '../../modules/igv/igv_report1'




workflow {
	
	if( params.help ) {
		log.info(paramsHelp())
		exit 0
		}  else {
		// Print summary of supplied parameters
		log.info paramsSummaryLog(workflow)
		}

	
		
	versions = Channel.empty()
	multiqc = Channel.empty()
	
	if(params.fasta==null) {
      log.error("undefined --fasta");
	  exit -1
      }
	
	if(params.samplesheet==null) {
      log.error("undefined --samplesheet");
	  exit -1
      }

	def metadata= [id:"igvreport"]


  /***************************************************
   *
   *  PREPARE FASTA REFERENCE
   *
   */
	PREPARE_ONE_REFERENCE(
			metadata.plus(scatter_bed:false),
			Channel.of(params.fasta).map{file(it)}.map{[[id:it.baseName],it]}
			)
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)

   /*READ BAM samplesheet */
    READ_BAMS(
        metadata.plus(arg_name:"samplesheet"),
        params.samplesheet
        )
	versions = versions.mix(READ_BAMS.out.versions)

	META_TO_BAMS(
			metadata ,
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
			READ_BAMS.out.samplesheet
			)
	versions = versions.mix(META_TO_BAMS.out.versions)


	if(params.variants!=null) {
		READ_VARIANTS(
			metadata.plus(arg_name:"variants"),
        	params.variants
			)
		versions = versions.mix(READ_VARIANTS.out.versions)

		ch1 = READ_VARIANTS.out.samplesheet
			.map{assertKeyExistsAndNotEmpty(it,"contig")}
			.map{assertKeyExistsAndNotEmpty(it,"sample")}
			.map{
				if(it.start!=null && it.end!=null) return it;
				if(it.position==null) throw new IllegalArgumentException("missing 'position' in ${meta}");
				def pos = (it.position as int)
				return it.plus([start:pos,end:pos]);
				}
			/** multiple sample name specified  with comma */
			.flatMap{row->
				def samples = row.sample.split(",");
				def L = [];
				def i;
				for(i=0;i< samples.size();i++) {
					def sn = samples[i].trim();
					if(sn.isEmpty()) continue;
					L.add(row.plus(sample:sn, sample_idx:i));
					}
				return L;
				}
			.combine(META_TO_BAMS.out.bams)
			.filter{meta1,meta2,bam,bai->meta1.sample==meta2.id || meta1.sample=="*"}
			.map{meta1,meta2,bam,bai->[[contig:meta1.contig,start:(meta1.start as int),end:(meta1.end as int)],[meta1,bam,bai]]}
			.groupTuple()
			// keep bam order as defined by 'sample_idx'
			.map{key_pos,arrays->[key_pos,arrays.sort{A1,A2->A1[0].sample_idx <=> A2[0].sample_idx } ]}
			.map{key_pos,arrays->[
				arrays[0][0].plus(id: key_pos.contig+"_"+key_pos.start+(key_pos.start==key_pos.end?"":"_"+key_pos.end)) /* first meta */]
				arrays.collect{A->A[1]}, //bams
				arrays.collect{A->A[2]} // bais
				}
		
		DOWNLOAD_CYTOBAND(PREPARE_ONE_REFERENCE.out.dict)
		versions = versions.mix(DOWNLOAD_CYTOBAND.out.versions)


		DOWNLOAD_REFGENE(PREPARE_ONE_REFERENCE.out.dict)
		versions = versions.mix(DOWNLOAD_REFGENE.out.versions)



		IGV_REPORT(
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
			PREPARE_ONE_REFERENCE.out.dict,
			DOWNLOAD_CYTOBAND.out.bed,
			DOWNLOAD_REFGENE.out.tabix,
			[id:"nosamplelist",[]],
			[id:"noped",[]],
			vcf_ch,
			ch1
			)
		versions = versions.mix(IGV_REPORT.out.versions)
		}




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
