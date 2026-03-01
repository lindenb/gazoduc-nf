<?xml version='1.0'  encoding="ISO-8859-1"?>
<xsl:stylesheet xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
	xmlns:date="http://exslt.org/dates-and-times" 
	version='1.0' 
	>
<xsl:output method="text"/>

<xsl:template match="/">
/*

THIS FILE WAS GENERATED DO NOT EDIT !!!

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

<xsl:apply-templates select="samplesheet"/>
</xsl:template>


<xsl:template match="samplesheet">

include { runOnComplete;dumpParams                 } from '../../modules/utils/functions.nf'
include { isBlank                                  } from '../../modules/utils/functions.nf'
include { parseBoolean                             } from '../../modules/utils/functions.nf'
include { PREPARE_ONE_REFERENCE                    } from '../../subworkflows/samtools/prepare.one.ref'
include { ENCODE_BLACKLIST                         } from '../../modules/encode/blacklist' 
include { BEDTOOLS_SUBTRACT                        } from '../../modules/bedtools/subtract' 
include { BEDTOOLS_INTERSECT                       } from '../../modules/bedtools/intersect' 
include { BED_CLUSTER                              } from '../../modules/jvarkit/bedcluster'
include { HAPLOTYPECALLER                          } from '../../subworkflows/gatk/haplotypecaller'
include { BEDTOOLS_MAKEWINDOWS                     } from '../../modules/bedtools/makewindows'
include { JVARKIT_VCFFILTERJDK                     } from '../../modules/jvarkit/vcffilterjdk'
include { BCFTOOLS_INDEX                           } from '../../modules/bcftools/index'
include { MERGE_VCFS                               } from '../../modules/gatk/mergevcfs'


<xsl:for-each select="bam[@status='control']">include{ VS_CONTROL as  <xsl:apply-templates select="." mode="vs_control"/> } from './sub.ultrarares.nf' 
</xsl:for-each>


workflow {
  versions = Channel.empty()
  multiqc = Channel.empty()


	def metadata = [id:"ultrares"]
	versions = Channel.empty()
	multiqc  = Channel.empty()


	if(params.fasta==null) {
			throw new IllegalArgumentException("undefined --fasta");
			}
	
    if(params.jvarkit_filter==null) {
			throw new IllegalArgumentException("undefined --jvarkit_filter");
			}
  
  /***************************************************
   *
   *  PREPARE FASTA REFERENCE
   *
   */
	PREPARE_ONE_REFERENCE(
			metadata,
			Channel.of(params.fasta).map{file(it)}.map{[[id:it.baseName],it]}
			)
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)

	if(params.vcf==null) {

		/***************************************************
		*
		*  MAKE BED
		*
		*/
		if(params.bed==null) {
			bed = PREPARE_ONE_REFERENCE.out.scatter_bed
			}
		else
			{
			BEDTOOLS_INTERSECT(
				PREPARE_ONE_REFERENCE.out.fai,
				Channel.of(params.bed)
					.map{bed->file(bed)}
					.map{bed->[[id:bed.baseName],bed]}
					.combine(PREPARE_ONE_REFERENCE.out.scatter_bed)
					.map{meta1,bed1,meta2,bed2->[meta1,bed1,bed2]}
				)
			versions = versions.mix(BEDTOOLS_INTERSECT.out.versions)
			bed = BEDTOOLS_INTERSECT.out.bed
			}


		/***************************************************
		*
		* Download encode black list
		*
		*/
		ENCODE_BLACKLIST(PREPARE_ONE_REFERENCE.out.dict)
		versions = versions.mix(ENCODE_BLACKLIST.out.versions)
		BEDTOOLS_SUBTRACT(bed.combine(ENCODE_BLACKLIST.out.bed).map{meta1,bed1,_meta2,bed2->[meta1,bed1,bed2]})
		versions = versions.mix(BEDTOOLS_SUBTRACT.out.versions)
		bed = BEDTOOLS_SUBTRACT.out.bed.first()


		BEDTOOLS_MAKEWINDOWS( bed )
	versions = versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions)


		/* if it's an exome , group the small genome together in BED */
		BED_CLUSTER(
			PREPARE_ONE_REFERENCE.out.dict,
			BEDTOOLS_MAKEWINDOWS.out.bed
			)
		versions = versions.mix(BED_CLUSTER.out.versions)
		bed = BED_CLUSTER.out.bed
			.map{_meta,bed->bed}
			.map{it instanceof List?it:[it]}
			.flatMap()
			.map{bed->[[id:bed.baseName],bed]}
		

		/***************************************************
		*
		* Run Hapcaller on CASES
		*
		*/
	HAPLOTYPECALLER(
		metadata.plus(
			gvcf_merge_method : "combinegvcfs",
			with_split_bed: false
			),
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		[[id:"oneref"],[]],
		[[id:"nodbsnp"],[],[]],
		[[id:"noped"],[]],
		bed,
		Channel.of(
		<xsl:for-each select="bam[@status='case']">
			<xsl:if test="position()&gt;1"><xsl:text>,
			</xsl:text>
			</xsl:if>
			<xsl:apply-templates select="." mode="channel"/>
		</xsl:for-each>)
		)
	versions = versions.mix(HAPLOTYPECALLER.out.versions)
	multiqc = multiqc.mix(HAPLOTYPECALLER.out.multiqc)
	vcf = 	HAPLOTYPECALLER.out.vcf
	}
else
	{
	vcf = Channel.of([[id:"user_vcf"],file(params.vcf),file(params.vcf+".tbi")])
	}
  

  jvarkit_filter = [[id:"jvarkit_filter"],file(params.jvarkit_filter)]
	JVARKIT_VCFFILTERJDK(
		jvarkit_filter,
		[[id:"no_ped"],[]],
		vcf
		)
	versions = versions.mix(JVARKIT_VCFFILTERJDK.out.versions)
	
    BCFTOOLS_INDEX(JVARKIT_VCFFILTERJDK.out.vcf)
    versions = versions.mix(BCFTOOLS_INDEX.out.versions)

	vcf = BCFTOOLS_INDEX.out.vcf

<xsl:for-each select="bam[@status='control']"> 
	<xsl:apply-templates select="." mode="vs_control"/>(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
    	PREPARE_ONE_REFERENCE.out.fai,
    	PREPARE_ONE_REFERENCE.out.dict,
		jvarkit_filter,
		vcf,
		Channel.of(<xsl:apply-templates select="." mode="channel"/>)
		)
	vcf = <xsl:apply-templates select="." mode="vs_control"/>.out.vcf
	versions = versions.mix(<xsl:apply-templates select="." mode="vs_control"/>.out.versions)
	multiqc = multiqc.mix(<xsl:apply-templates select="." mode="vs_control"/>.out.multiqc)
	
</xsl:for-each>

  MERGE_VCFS(
	vcf.map{meta,vcf,tbi->[vcf,tbi]}
		.flatMap()
		.collect()
		.map{files->[metadata,files.sort()]}
  	)
  versions = versions.mix(MERGE_VCFS.out.versions)
}
</xsl:template>

<xsl:template match="bam" mode="channel">
  <xsl:text>[[id:"</xsl:text>
  <xsl:value-of select="@sample"/>
  <xsl:text>",status:"</xsl:text>
  <xsl:value-of select="@status"/>
  <xsl:text>"], file("</xsl:text>
  <xsl:value-of select="@path"/>
  <xsl:text>") ,file("</xsl:text>
  
  <xsl:choose>
	<xsl:when test=" substring(@path, string-length(@path) - 3) = '.bam'">
		<xsl:value-of select="@path"/>
		<xsl:text>.bai</xsl:text>
		</xsl:when>
	<xsl:when test=" substring(@path, string-length(@path) - 4) = '.cram'">
		<xsl:value-of select="@path"/>
		<xsl:text>.crai</xsl:text>
	</xsl:when>
	<xsl:otherwise>
		<xsl:text>boum_bam_suffix: </xsl:text>
		<xsl:value-of select="substring(@path, string-length(@path) - 5)"/>
	</xsl:otherwise>
  </xsl:choose>
  <xsl:text>")]</xsl:text>
</xsl:template>

<xsl:template match="bam" mode="vs_control">
<xsl:text>VS_CONTROL_</xsl:text>
 <xsl:value-of select="translate(@sample,'-.', '_')"/>
</xsl:template>


</xsl:stylesheet>

