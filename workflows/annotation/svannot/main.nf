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


include { validateParameters                          } from 'plugin/nf-schema'
include { paramsHelp                                  } from 'plugin/nf-schema'
include { paramsSummaryLog                            } from 'plugin/nf-schema'
include { samplesheetToList                           } from 'plugin/nf-schema'
include { runOnComplete;dumpParams                    } from '../../../modules/utils/functions.nf'
include { PREPARE_ONE_REFERENCE                       } from '../../../subworkflows/samtools/prepare.one.ref'
include { VCF_INPUT                                   } from '../../../subworkflows/nf/vcf_input'
include { DOWNLOAD_GTF_OR_GFF3  as  DOWNLOAD_GTF      } from '../../../modules/gtf/download'
include { DOWNLOAD_GNOMAD_SV                          } from '../../../modules/gnomad_sv/download.vcf'
include { GTF_TO_BED                                  } from '../../../modules/jvarkit/gtf2bed'
include { JVARKIT_VCFGNOMADSV                         } from '../../../modules/jvarkit/vcfgnomadsv'
include { BEDTOOLS_SLOP                               } from '../../../modules/bedtools/slop'
include { BEDTOOLS_MERGE                              } from '../../../modules/bedtools/merge'
include { BCFTOOLS_VIEW                               } from '../../../modules/bcftools/view'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_INV          } from '../../../modules/bcftools/view'
include { ANNOTSV                                     } from '../../../subworkflows/annotsv'
include { DGV_DOWNLOAD                                } from '../../../modules/dgv/download'
include { READ_SAMPLESHEET                            } from '../../../subworkflows/nf/read_samplesheet'
include { META_TO_BAMS                                } from '../../../subworkflows/samtools/meta2bams1'
include { BCFTOOLS_QUERY                              } from '../../../modules/bcftools/query'
include { COVERAGE_GRID                               } from '../../../modules/jvarkit/coveragegrid'
include { GHOSTSCRIPT_MERGE                           } from '../../../modules/gs/merge/main.nf'
include { PDF_NAVIGATION                              } from '../../../modules/pdf/navigation/main.nf'
workflow {
		validateParameters()

		if( params.help ) {
			log.info(paramsHelp())
			exit 0
		}  else {
		// Print summary of supplied parameters
		log.info paramsSummaryLog(workflow)
		}



	if(params.fasta==null) {
			throw new IllegalArgumentException("undefined --fasta");
			}
	if(params.samplesheet==null) {
		throw new IllegalArgumentException("undefined --samplesheet");
		}
	multiqc = Channel.empty()
	versions = Channel.empty()
	def metadata = [id:"annotsv"]
		

	PREPARE_ONE_REFERENCE(
		metadata,
		Channel.fromPath(params.fasta).map{f->[[id:f.baseName],f]}
		)
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)
	

	READ_SAMPLESHEET(
        [arg_name:"samplesheet"],
        params.samplesheet
        )
	versions = versions.mix(READ_SAMPLESHEET.out.versions)

	META_TO_BAMS(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		READ_SAMPLESHEET.out.samplesheet
		)
	versions = versions.mix(META_TO_BAMS.out.versions)




	DOWNLOAD_GTF(PREPARE_ONE_REFERENCE.out.dict)
	versions = versions.mix(DOWNLOAD_GTF.out.versions)

    DOWNLOAD_GNOMAD_SV( PREPARE_ONE_REFERENCE.out.dict )
    versions = versions.mix(DOWNLOAD_GNOMAD_SV.out.versions)

    DGV_DOWNLOAD( PREPARE_ONE_REFERENCE.out.dict )
    versions = versions.mix(DGV_DOWNLOAD.out.versions)


	VCF_INPUT(
        metadata.plus([
			path:params.vcf,
			require_index:true,
			arg_name:"vcf",
			required: true,
			unique : false
			])
        )
    versions = versions.mix(VCF_INPUT.out.versions)

	if(params.bed==null) {
		bed = PREPARE_ONE_REFERENCE.out.scatter_bed
	} else if(params.bed.matches("(coding[_-])?(exon|gene)s?")) {
		GTF_TO_BED(
			PREPARE_ONE_REFERENCE.out.dict,
			DOWNLOAD_GTF.out.gtf.map{meta,gtf,_tbi->[meta,gtf]}
			)
		versions = versions.mix(GTF_TO_BED.out.versions)
		bed =  GTF_TO_BED.out.bed
	} else {
		bed =Channel.of(params.bed).map{[[id:"bed"],file(f)]}
	}

	BEDTOOLS_SLOP(
		PREPARE_ONE_REFERENCE.out.fai,
		bed
		)
	versions = versions.mix(BEDTOOLS_SLOP.out.versions)

	BEDTOOLS_MERGE( BEDTOOLS_SLOP.out.bed )
	versions = versions.mix(BEDTOOLS_MERGE.out.versions)

	BCFTOOLS_VIEW(
		BEDTOOLS_MERGE.out.bed,
		[[id:"no_sample"],[]],
		VCF_INPUT.out.vcf
		)
	versions = versions.mix(BCFTOOLS_VIEW.out.versions)

	JVARKIT_VCFGNOMADSV(
		DOWNLOAD_GNOMAD_SV.out.vcf,
		BCFTOOLS_VIEW.out.vcf
		)
	versions = versions.mix(JVARKIT_VCFGNOMADSV.out.versions)

	ANNOT_WITH_GTF(
		PREPARE_ONE_REFERENCE.out.fai,
		DOWNLOAD_GTF.out.gtf,
		JVARKIT_VCFGNOMADSV.out.vcf
		)
	versions = versions.mix(ANNOT_WITH_GTF.out.versions)


	ANNOTSV(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		JVARKIT_VCFGNOMADSV.out.vcf
		)
	versions = versions.mix(ANNOTSV.out.versions)


	BCFTOOLS_VIEW_INV(
		[[id:"nobed"],[]],
		[[id:"no_sample"],[]],
		JVARKIT_VCFGNOMADSV.out.vcf.map{meta,vcf->[meta,vcf,[]]}
		)
	versions = versions.mix(BCFTOOLS_VIEW_INV.out.versions)


	MAKE_SAMPLESHEET2(
		META_TO_BAMS.out.bams
			.map{meta,bam,bai->"${meta.id}\tBAMS/${bam.name}\t${meta.status=="case"?"red":"green"}"}
			.collect()
			.map{L->[[id:"bams"],L.sort()]}
	)
	
	BCFTOOLS_QUERY(BCFTOOLS_VIEW_INV.out.vcf)
	versions = versions.mix(BCFTOOLS_QUERY.out.versions)
	
	

	COVERAGE_GRID(
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		DOWNLOAD_GTF.out.gtf,
		DOWNLOAD_GNOMAD_SV.out.vcf,
		MAKE_SAMPLESHEET2.out.samplesheet,
		META_TO_BAMS.out.bams
			.map{meta,bam,bai->[bam,bai]}
			.flatMap()
			.collect()
			.map{L->[[id:"bams"],L.sort()]},
		BCFTOOLS_QUERY.out.output
			.map{_meta,f->f}
			.splitCsv(header:false,sep:'\t')
			.map{row->[[id:"${row[0]}_${row[1]}_${row[2]}",title:"${row[0]}:${row[1]}-${row[2]}"],row[0],(row[1] as int),(row[2] as int)]}
		)
	versions = versions.mix(COVERAGE_GRID.out.versions)

	GHOSTSCRIPT_MERGE(COVERAGE_GRID.out.postscript)
	versions = versions.mix(GHOSTSCRIPT_MERGE.out.versions)

	PDF_NAVIGATION(GHOSTSCRIPT_MERGE.out.pdf.map{m,pdf->pdf}.collect().map{L->[[id:"inv"],L.sort()]})
	versions = versions.mix(PDF_NAVIGATION.out.versions)
}

process ANNOT_WITH_GTF {
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(fai)
	tuple val(meta2),path(gtf),path(gtf_tbi)
	tuple val(meta3),path(genes_of_interest)
	tuple val(meta ),path(vcf)
output:
	tuple val(meta),path("*.vcf.gz"),emit:vcf
script:
	def xstream = task.ext.xstream?:"1000"
"""
mkdir -p TMP

gunzip -c "${gtf}" |\\
	awk -F '\t' '(\$3=="exon") {printf("%s\t%d\t%s\\n",\$1,int(\$2)-1,int(\$3))}' |\\
	sort -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n -T TMP |\\
	bedtools merge |\\
	sed 's/\$/\t1/' > TMP/exons.bed

gunzip -c "${gtf}" |\\
	awk -F '\t' '(\$3=="CDS") {printf("%s\t%d\t%s\\n",\$1,int(\$2)-1,int(\$3))}' |\\
	sort -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n -T TMP |\\
	bedtools merge |\\
	sed 's/\$/\t1/' > TMP/cds.bed

gunzip -c "${gtf}" |\\
	awk -F '\t' '(\$3=="gene") {printf("%s\t%d\t%s\\n",\$1,int(\$2)-1,int(\$3))}' |\\
	sort -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n -T TMP |\\
	bedtools merge |\\
	sed 's/\$/\t1/' > TMP/genes.bed

gunzip -c "${gtf}" |\\
	awk -F '\t' '(\$3=="transcript") {printf("%s\t%d\t%s\\n",\$1,int(\$2)-1,int(\$3))}' |\\
	sort -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n -T TMP |\\
	bedtools merge |\\
	sed 's/\$/\t1/' > TMP/transcript.bed

gunzip -c "${gtf}" |\\
	grep -F -w protein_coding |\\
	awk -F '\t' '(\$3=="gene") {printf("%s\t%d\t%s\\n",\$1,int(\$2)-1,int(\$3))}' |\\
	sort -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n -T TMP |\\
	bedtools merge |\\
	sed 's/\$/\t1/' > TMP/protein_coding.bed

bedtools subtract \\
	-a TMP/transcript.bed \\
	-b TMP/exons.bed |\\
	cut -f1,2,3 |\\
	sort -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n -T TMP |\\
	bedtools merge |\\
	sed 's/\$/\t1/' > TMP/introns.bed

gunzip -c "${gtf}" |\\
	awk -F '\t' '(\$3=="gene") {S=\$7; B=int(\$2)-1;E=int(\$3);if(S=="+") {E=B-1 ; B = B-${xstream}; if(B<0)B=0; if(E<0) E=0;} else {B=E;E=E+${xstream};} printf("%s\t%d\t%s\\n",\$1,B,E)}' |\\
	sort -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n -T TMP |\\
	bedtools merge |\\
	sed 's/\$/\t1/' > TMP/upstream.bed

gunzip -c "${gtf}" |\\
	awk -F '\t' '(\$3=="gene") {S=\$7; B=int(\$2)-1;E=int(\$3);if(S=="-") {E=B-1 ; B = B-${xstream}; if(B<0)B=0; if(E<0) E=0;} else {B=E;E=E+${xstream};} printf("%s\t%d\t%s\\n",\$1,B,E)}' |\\
	sort -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n -T TMP |\\
	bedtools merge |\\
	sed 's/\$/\t1/' > TMP/downstream.bed

if ${genes_of_interest?true:false}
then

	cat "${genes_of_interest}" |\\
		grep -vE '^#' |\\
		grep -vE '^\$' |\\
		sort -S ${task.memory.kilo} -t '\t' -k1,1 -T TMP > TMP/jeter.a

	gunzip -c "${gtf}" |\\
		awk -F '\t' '(\$3=="gene") |\\
		java -jar \${HOME}/jvarkit.jar gtf2bed -c 'gene_name' |\\
		sort -S ${task.memory.kilo} -t '\t' -k4,4 -T TMP > TMP/jeter.b

	join -t '\t' -1 1 -2 4 -o '2.1,2.2,2.3,2.4' | uniq > TMP/roi.bed

fi



"""
}

process MAKE_SAMPLESHEET2 {
	executor "local"
	input:
		tuple val(meta),val(L)
	output:
		tuple val(meta),path("samplesheet.tsv"),emit:samplesheet
	script:
"""
cat << EOF > samplesheet.tsv
sample\tbam\tcolor
${L.join("\n")}
EOF
"""
}