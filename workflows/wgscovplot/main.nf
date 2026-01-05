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
include { runOnComplete;dumpParams                 } from '../../modules/utils/functions.nf'
include { median                                   } from '../../modules/utils/functions.nf'
include { PREPARE_ONE_REFERENCE                    } from '../../subworkflows/samtools/prepare.one.ref'
include { META_TO_BAMS                             } from '../../subworkflows/samtools/meta2bams1'
include { READ_SAMPLESHEET                         } from '../../subworkflows/nf/read_samplesheet'
include { FAI_TO_SVG                               } from '../../subworkflows/circular/fai2svg'
include { CYTOBAND_TO_SVG                          } from '../../subworkflows/circular/cytoband2svg'
include { DICT_TO_BED                              } from '../../modules/jvarkit/dict2bed'
include { BEDTOOLS_MAKEWINDOWS                     } from '../../modules/bedtools/makewindows'
include { BEDTOOLS_INTERSECT                       } from '../../modules/bedtools/intersect'
include { BED_TO_XML                               } from '../../modules/jvarkit/bed2xml'
include { MOSDEPTH                                 } from '../../modules/mosdepth'
include { XSLTPROC                                 } from '../../modules/xsltproc'
include { BATIK                                    } from '../../subworkflows/batik'

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
		metadata.plus(radius: (params.radius as double)),
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

		CYTOBAND_TO_SVG(
			metadata,
			FAI_TO_SVG.out.dict,
			Channel.empty()
			)
		versions = versions.mix(CYTOBAND_TO_SVG.out.versions)

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
			
			
		max_cov = mean_cov_ch
			.map{meta,cov->cov}
			.collect()
			.map{L->L.max()}
			

		regions_ch = MOSDEPTH.out.regions_bed.join(mean_cov_ch)
			.map{meta,bed,_csi,cov->[meta.plus(coverage:cov),bed]}
			.combine(max_cov)
			.map{meta,bed,max_cov->[meta.plus(coverage_max:max_cov),bed]}
			.collect(flat:false)
			.flatMap{rows->{
				def L=rows.collect();
				def r_max = (metadata.radius*1.0)-50;
				def r_min  = r_max*0.1;
				def dr = java.lang.Math.min(20.0,(r_max-r_min)/L.size());
				def r = r_max;
				L = L.sort{A,B-> A[0].coverage <=>  B[0].coverage}
				def i=0;
				while( i < L.size() ) {
					def item = L[i];
					item[0] = item[0].plus(
						index : i,
						index_mod2 : (i%2),/* cannot make index%2 work in config... */
						count: L.size(),
						outer_radius: r,
						inner_radius: r - dr *0.95
						) 
					r -= dr;
					i++;
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


		XSLTPROC(
			[[id:"xslt"],file("${moduleDir}/../../src/xsl/histogram.xsl")],
			BED_TO_XML.out.xml
			)
		versions = versions.mix(XSLTPROC.out.versions)

		BUILD_SVG(
			MOSDEPTH_TO_BED.out.css
				.map{meta,css->css}
				.collect()
				.map{files->[[id:"css"],files.sort()]},
			MOSDEPTH_TO_BED.out.defs
				.map{meta,svg->svg}
				.collect()
				.map{files->[[id:"defs"],files.sort()]},
			FAI_TO_SVG.out.svg,
			CYTOBAND_TO_SVG.out.svg,
			XSLTPROC.out.xml
				.map{meta,svg->svg}
				.collect()
				.map{files->[[id:"svgs"],files.sort()]}
			)
		versions = versions.mix(BUILD_SVG.out.versions)

		BATIK(
			metadata,
			BUILD_SVG.out.svg
			)
		versions = versions.mix(BATIK.out.versions)
	emit:
		versions
		multiqc
		svg  = BUILD_SVG.out.svg
	}


process MOSDEPTH_TO_BED {
label "process_single"
tag "${meta.id}"

input:
	tuple val(meta),path(bed)
output:
	tuple val(meta),path("*.bed"),emit:bed
	tuple val(meta),path("*.css"),emit:css
	tuple val(meta),path("*.defs.svg"),emit:defs
	path("versions.yml"),emit:versions
script:
	def cov = meta.coverage as double
	def max_cov = meta.coverage_max as double
	def prefix = task.ext.prefix?:"${meta.id}"
	def md5 = "${meta.id.md5().substring(0,5)}.histo"
	def outer_radius = meta.outer_radius as double
	def inner_radius = meta.inner_radius as double
	def index = meta.index?:0
"""
${bed.name.endsWith(".gz")?"gunzip -c":"cat"} ${bed} |\\
	awk -F '\t' '{
		COV = \$4 * 1.0;
		printf("%s\t%s\t%s\t%f",\$1,\$2,\$3,COV);
		printf("\t%f",  ( COV / ${max_cov}) );
		if( COV >= ${cov} * 2.0) {
			printf("\tbar DUP");
			}
		else if( COV <= ${cov} * 0.5) {
			printf("\tbar DEL");
			}
		else
			{
			printf("\tbar PASS");
			}
		printf("\\n");
		}' > ${prefix}.bed

# gradient doesn't work, should use a mask for each bar.. borrrrring
touch ${prefix}.css
touch ${prefix}.defs.svg

cat << EOF > versions.yml
${task.process}:
	awk: \$(awk --version | awk '(NR==1)')
EOF
"""
stub:
	def prefix = task.ext.prefix?:"${meta.id}"
"""
touch versions.yml ${prefix}.bed ${prefix}.css ${prefix}.defs.svg
"""
}


process BUILD_SVG {
tag "${meta.id}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(css)
    tuple val(meta2),path(defs)
	tuple val(meta3),path(fai_svg)
	tuple val(meta4),path(cytoband_svg)
	tuple val(meta ),path(fragments)
output:
    tuple val(meta),path("*.svg"),emit:svg
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"${meta.id}.circular"
    def radius = (task.ext.radius?:1000) as int
    def xmlns = "http://www.w3.org/2000/svg"
    def margin = (task.ext.margin?:100) as int
    def margin_top = (task.ext.margin_top?:margin) as int
    def margin_left = (task.ext.margin_left?:margin) as int
    def margin_bottom = (task.ext.margin_bottom?:margin) as int
    def margin_right = (task.ext.margin_right?:margin) as int
    def diameter = (radius*2)
    def image_width = diameter + margin_left + margin_right
    def image_height = diameter + margin_top + margin_bottom
    def title = task.ext.title?:"${meta.title?:meta.id}"
"""
touch versions.yml
mkdir -p TMP

cat << EOF > jeter.css

rect.bckg {
	stroke:none;
	fill:rgb(254,254,254);
	}

path.donut_0 {
		opacity: 0.8;
        stroke: gray;
        stroke-dasharray: 2, 2;
        fill: beige;
        }

path.donut_1 {
		opacity: 0.8;
        stroke: gray;
        stroke-dasharray: 2, 2;
        fill: blanchedalmond;
        }

path.bar {
	opacity: 0.8;
	 stroke:gray;
	stroke-width:0.2;
	fill:lightgray;
	}

path.DEL {
	fill:red;
	}
path.DUP {
	fill:blue;
	}

line.cross {
	stroke:gray;
	stroke-width:1;
	}

EOF
cat ${css} >> jeter.css
cat "${moduleDir}/../../src/css/fai.css" >> jeter.css

echo '<defs xmlns="${xmlns}">' > jeter.defs.xml
cat ${defs}  >> jeter.defs.xml
echo '</defs>'  >> jeter.defs.xml


echo '<g xmlns="${xmlns}">' > jeter.fragments.xml
cat ${fragments}  >> jeter.fragments.xml
echo '</g>'  >> jeter.fragments.xml

cat << __EOF__ > jeter.doc.xml
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "https://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg
    version="1.1" 
    xmlns="${xmlns}"
    xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
    xmlns:foaf="http://xmlns.com/foaf/0.1/"
    xmlns:dc="http://purl.org/dc/elements/1.1/"
    xmlns:xi="http://www.w3.org/2001/XInclude"
    width="${image_width +1}" height="${image_height +1}">
  <metadata>
    <rdf:RDF>
        <rdf:Description rdf:about="">
            <dc:author>Pierre Lindenbaum</dc:author>
        </rdf:Description>
    </rdf:RDF>
  </metadata>
  <style>
  <xi:include parse="text" href="jeter.css"/>
  </style>
  <xi:include href="jeter.defs.xml"/>
  <g>
    <rect x="0" y="0" width="${image_width}" height="${image_height}" class="bckg"/>
    <g transform="translate(${margin_left},${margin_top})">
        <g transform="translate(${radius},${radius})">
			<xi:include href="${fai_svg.name}"/>
			<xi:include href="${cytoband_svg.name}"/>
            <xi:include href="jeter.fragments.xml"/>
            <!-- central cross -->
            <text x="0" y="0" class="title"><![CDATA[${title}]]></text>
            <g>
                <line x1="-10" y1="0" x2="10"  y2="0" class="cross"/>
                <line y1="-10" x1="0" y2="10"  x2="0" class="cross"/>
            </g>
        </g>
    </g>
   
  </g>
</svg>
__EOF__


xmllint --xinclude jeter.doc.xml > TMP/jeter2.svg
mv TMP/jeter2.svg TMP/jeter.svg


mv TMP/jeter.svg "${prefix}.svg"

### gzip --best < "${prefix}.svg" > \${HOME}/jeter.svg.gz
"""
stub:
    def prefix = task.ext.prefix?:"${meta.id}.circular"
"""
touch versions.yml  ${prefix}.svg
"""
}

