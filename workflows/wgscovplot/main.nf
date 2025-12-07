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
nextflow.enable.dsl=2
include { runOnComplete;dumpParams                 } from '../../modules/utils/functions.nf'
include { median                                   } from '../../modules/utils/functions.nf'
include { PREPARE_ONE_REFERENCE                    } from '../../subworkflows/samtools/prepare.one.ref'
include { META_TO_BAMS                             } from '../../subworkflows/samtools/meta2bams1'
include { READ_SAMPLESHEET                         } from '../../subworkflows/nf/read_samplesheet'
include { FAI_TO_SVG                               } from '../../subworkflows/circular/fai2svg'
include { DICT_TO_BED                              } from '../../modules/jvarkit/dict2bed'
include { BEDTOOLS_MAKEWINDOWS                     } from '../../modules/bedtools/makewindows'
include { BEDTOOLS_INTERSECT                       } from '../../modules/bedtools/intersect'
include { BED_TO_XML                               } from '../../modules/jvarkit/bed2xml'
include { MOSDEPTH                                 } from '../../modules/mosdepth'

workflow {

	def metadata = [id:"ultrares"]
	versions = Channel.empty()
	multiqc  = Channel.empty()


	if(params.fasta==null) {
			throw new IllegalArgumentException("undefined --fasta");
			}
	if(params.samplesheet==null) {
			throw new IllegalArgumentException("undefined --samplesheet");
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


	PLOT_WGS_COVERAGE(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		PREPARE_ONE_REFERENCE.out.scatter_bed,
		META_TO_BAMS.out.bams
		)
	versions = versions.mix(PLOT_WGS_COVERAGE.out.versions)
	}

runOnComplete(workflow)

workflow PLOT_WGS_COVERAGE {
	take:
		metadata
		fasta
		fai
		dict
		scatter_bed
		bams
	main:
		versions = Channel.empty()
		multiqc  = Channel.empty()

		FAI_TO_SVG(metadata,fai,dict)
		versions = versions.mix(FAI_TO_SVG.out.versions)

		DICT_TO_BED(FAI_TO_SVG.out.dict)
		versions = versions.mix(DICT_TO_BED.out.versions)

		BEDTOOLS_INTERSECT(
			fai,
			DICT_TO_BED.out.bed.combine(scatter_bed)
				.map{meta1,bed1,_meta2,bed2->[meta1,bed1,bed2]}
			)
		versions = versions.mix(BEDTOOLS_INTERSECT.out.versions)


		genome_length_ch = BEDTOOLS_INTERSECT.out.bed
            .splitCsv(sep:'\t',header:false)
            .map{meta,tokens->{
                return [meta,(tokens[2] as long)-(tokens[1] as long)];
                }
                }
            .groupTuple()
            .map{meta,lengths->[meta,lengths.sum()]}

		bed_ch =  BEDTOOLS_INTERSECT.out.bed
			.join(genome_length_ch)
			.map{meta,bed,genome_size->[meta.plus(genome_length:genome_size),bed]}

		BEDTOOLS_MAKEWINDOWS(bed_ch)
		versions = versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions)

		MOSDEPTH(
			fasta,
			fai,
			bams.combine(BEDTOOLS_MAKEWINDOWS.out.bed)
					.map{meta1,bam,bai,_meta2,bed->[meta1,bam,bai,bed]}
			)
		versions = versions.mix(MOSDEPTH.out.versions)
		mean_cov_ch = MOSDEPTH.out.summary_txt
			.splitCsv(sep:'\t',header:true)
			.filter{_meta,row->row.chrom=="total_region"}
            .map{meta,row->[meta,java.lang.Math.max(1.0,(row.mean as double))]}
            .groupTuple()
			.map{meta,array->[meta,array[0]]}
			.view()
			
		max_cov = mean_cov_ch
			.map{meta,cov->cov}
			.collect()
			.map{L->L.max()}
			.view()

		regions_ch = MOSDEPTH.out.regions_bed.join(mean_cov_ch)
			.map{meta,bed,csi,cov->[meta.plus(coverage:cov),bed,csi]}
			.view()
			.combine(max_cov)
			.map{meta,bed,csi,max_cov->[meta.plus(coverage_max:max_cov),bed,csi]}
			.collect(flat:false)
			.flatMap{rows->{
				def L=rows.collect();
				log.warn("rows== ${rows}");
				log.warn("L == ${L}");
				L = L.sort{A,B-> A[0].coverage <=>  B[0].coverage}
				def i=0;
				for(i=0;i< L.size();i++) {
					def item = L[i];
					item[0] = item[0].plus(
						index : i,
						count: L.size()
						) 
					}				
				return L;
				}}
		MOSDEPTH_TO_BED(regions_ch)
		versions = versions.mix(MOSDEPTH_TO_BED.out.versions)
		
		BED_TO_XML(
			FAI_TO_SVG.out.dict.first(),
			MOSDEPTH_TO_BED.out.bed
			)
		versions = versions.mix(BED_TO_XML.out.versions)
	emit:
		versions
		multiqc
	}


process MOSDEPTH_TO_BED {
input:
	tuple val(meta),path(bed)
output:
	tuple val(meta),path("*.bed"),emit:bed
	path("versions.yml"),emit:versions
script:
	def cov = meta.coverage as double
	def max_cov = meta.coverage_max as double
	def prefix = task.ext.prefix?:"${meta.id}"
"""
${bed.name.endsWith(".gz")?"gunzip -c":"cat"} ${bed} |\\
	awk -F '\t' '{
		NORM = ( \$4 / ${cov})* ${max_cov};
		printf("%s\t%s\t%s\t%s\t%f\\n",\$1,\$2,\$3, \$4, NORM);
		}' > ${prefix}.bed
touch versions.yml
"""
stub:
	def prefix = task.ext.prefix?:"${meta.id}"
"""
touch versions.yml ${prefix}.bed
"""
}


process ZIPIT {
tag "N=${L.size()}"
input:
	val(meta)
	val(L)
output:
	path("${params.prefix?:""}plots.zip"),emit:zip
script:
"""
hostname 1>&2

cat << __EOF__ >  "${params.prefix?:""}index.html" 
<html>
<head>
<meta charset="UTF-8"/>
<title>${params.prefix?:""}.IndexCov</title>
<script>
var index=0;
var svgs=[${L.findAll{T->T.toString().endsWith(".svg")}.collect{T->"\""+file(T).name+"\""}.join(",")}];

function change(dx) {
	index+=dx;
	index = index%svgs.length;
	var E = document.getElementById("theimg");
	E.setAttribute("src",svgs[index]);
	}
</script>
</head>
<body>
<button onclick="change(-1)">Prev</button>
<button onclick="change( 1)">Next</button>
<br/>
<img id="theimg" src="${file(L[0]).name}" width="1000" />
</body>
</html>
__EOF__

cat << __EOF__ >  "${params.prefix?:""}all.html"
<html>
<head>
<meta charset="UTF-8"/>
<title>${params.prefix?:""}.WGSPlotCov</title>
<script>
function init() {
	var svgs=[${L.findAll{T->T.toString().endsWith(".svg")}.collect{T->"\""+file(T).name+"\""}.join(",")}];
	var i;
	var main=document.getElementById("main");
	for(i in svgs) {
		var img = document.createElement("img");
		img.setAttribute("src",svgs[i]);
		img.setAttribute("width","1000");
		main.appendChild(img);
		var br = document.createElement("br");
		main.appendChild(br);
		}
	}

window.addEventListener('load', (event) => {init();});

</script>
</head>
<body>
<div id="main"></div>
</body>
</html>
__EOF__


zip -j -9  "${params.prefix?:""}plots.zip" ${L.join(" ")} \
	"${params.prefix?:""}all.html" \
	"${params.prefix?:""}index.html"
"""
}
