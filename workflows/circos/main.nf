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

include { dumpParams;runOnComplete      } from '../../modules/utils/functions.nf'
include { isBlank                       } from '../../modules/utils/functions.nf'
include { PREPARE_ONE_REFERENCE         } from '../../subworkflows/samtools/prepare.one.ref'
include { BEDTOOLS_MAKEWINDOWS          } from '../../modules/bedtools/makewindows'
include { assertKeyExistsAndNotEmpty    } from '../../modules/utils/functions.nf'
include { assertKeyMatchRegex           } from '../../modules/utils/functions.nf'
include { DOWNLOAD_CYTOBAND             } from '../../modules/ucsc/download.cytobands'
include { BEDTOOLS_NUC                  } from '../../modules/bedtools/nuc'
include { SAMTOOLS_BEDCOV               } from '../../modules/samtools/bedcov'
include { CIRCOS                        } from '../../subworkflows/circos'

if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}

workflow {
	versions = Channel.empty()
	multiqc = Channel.empty()

	if(params.json==null) {
		throw new IllegalArgumentException("undefined --json");
		}
	if(params.fasta==null) {
		throw new IllegalArgumentException("undefined --fasta");
		}t


	def workflow_meta = [
		id: "circos"
		]

	PREPARE_ONE_REFERENCE(
		workflow_meta,
		Channel.fromPath(params.fasta).map{[[id:it.baseName],it]}
		)
    versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)
	
	dispatch = Channel.fromPath(params.json)
		.splitJson()
		.collect()
		.flatMap{v->
			def L=[];
			for(int i=0;i< v.size();i++) {
				def vi = v[i];
				if(isBlank(vi.order)) vi= vi.plus(order:java.lang.String.format("%03d",i));
				if(isBlank(vi.weight)) vi= vi.plus(weight:1.0);
				L.add(vi);
				}
			return L;
			}
		.branch{row->
			blank: isBlank(row.tracktype)
			bamcov: row.tracktype.equals("bamcov")
			gtf: row.tracktype.equals("gtf")
			other: true
		}
	

	BEDTOOLS_MAKEWINDOWS(PREPARE_ONE_REFERENCE.out.bed)
	versions = versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions)


	dispatch.blank
		.mix(dispatch.other)
		.map{throw new IllegalArgumentException("unknow or undefined trackType in ${it}")}
	

	//download cytobands
	DOWNLOAD_CYTOBAND(
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict
		)
	versions = versions.mix(DOWNLOAD_CYTOBAND.out.versions)

	// get GC% using bedtools nuc
	BEDTOOLS_NUC(
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		BEDTOOLS_MAKEWINDOWS.out.bed
		)
	versions = versions.mix(BEDTOOLS_NUC.out.versions)

	circos_files = Channel.empty()

	SAMTOOLS_BEDCOV(
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		BEDTOOLS_MAKEWINDOWS.out.bed, 
		dispatch.bamcov
			.map{assertKeyExistsAndNotEmpty(it,"bam")}
			.map{assertKeyExistsAndNotEmpty(it,"sample")}
			.map{assertKeyMatchRegex(it,"bam","^\\S+\\.(bam|cram)\$")}
			.map{isBlank(it.bai)?it.plus(bai: it.bam+(it.bam.endsWith(".bam")?".bai":".crai")):it}
			.map{assertKeyMatchRegex(it,"bai","^\\S+\\.(bai|crai)\$")}
			.map{isBlank(it.id)?it.plus(id: it.sample):it}
			.map{[
				it.findAll{k,v->!k.matches("bam|bai")},
				file(it.bam),
				file(it.bai)
				]}
		)
    versions = versions.mix(SAMTOOLS_BEDCOV.out.versions)
	circos_files = circos_files.mix(SAMTOOLS_BEDCOV.out.bed)

	CIRCOS(
		workflow_meta,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		BEDTOOLS_MAKEWINDOWS.out.bed, 
		DOWNLOAD_CYTOBAND.out.bed,
		circos_files
		)
	versions = versions.mix(CIRCOS.out.versions)
    }


runOnComplete(workflow);

