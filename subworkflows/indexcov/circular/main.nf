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

include { BED_TO_XML            } from '../../../modules/jvarkit/bed2xml'
include { XSLTPROC              } from '../../../modules/xsltproc'
include { BATIK                 } from '../../../subworkflows/batik'
include { FAI_TO_SVG            } from '../../../subworkflows/circular/fai2svg'
include { CYTOBAND_TO_SVG       } from '../../../subworkflows/circular/cytoband2svg'
include { DICT_TO_BED           } from '../../../modules/jvarkit/dict2bed'


workflow INDEXCOV_TO_SVG {
	take:
		metadata
		fasta
		fai
		dict
		scatter_N_bed
		bed
	main:
		versions  = Channel.empty()
		
		ch1 = bed
			.map{meta,bed->[meta,bed,bed]}//duplicate bed
			.splitCsv(header:false,sep:'\t',limit:1)
			.flatMap{meta,header,bed->{
				def L=[];
				def r_max = ((metadata.radius?:1000)*1.0)-50;
				def r_min  = r_max*0.1;
				def dr = java.lang.Math.min(50.0,(r_max-r_min)/(header.size()-3));
				def r = r_max;

				for(int i=3;i< header.size();i++) {
					def hash=[
						id: header[i],
						index : i-3,
						index_mod2 : ((i-3)%2),
						count : header.size()-3,
						file_id : meta.id,
						outer_radius: r,
                        inner_radius: r - dr *0.95
						]
					L.add([hash,bed]);
					r -= dr;
					}
				return L;
				}}

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

		SAMPLE_TO_BED(scatter_N_bed, ch1)
		versions = versions.mix(SAMPLE_TO_BED.out.versions)

		BED_TO_XML(FAI_TO_SVG.out.dict.first(), SAMPLE_TO_BED.out.bed)
		versions = versions.mix(BED_TO_XML.out.versions)
		
		XSLTPROC(
			[[id:"indexcov2svg"],"${moduleDir}/../../../src/xsl/indexcov2svg.xsl"],
			BED_TO_XML.out.xml
			)
    	versions = versions.mix(XSLTPROC.out.versions)

		BUILD_SVG(
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
	}


process SAMPLE_TO_BED {
label "process_single"
tag "${meta.id}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(scatter_N_bed)
	tuple val(meta ),path(bed)
output:
	tuple val(meta ),path("*.bed"),emit:bed
	path("versions.yml"),emit:versions
script:
	def column = ((meta.index as int)+4)
	def treshold = (task.ext.treshold?:0.1)
	def merge_d = task.ext.merge_d?:170000
	def min_size = merge_d*4
	def prefix = task.ext.prefix?:"${meta.id}.extract"
"""
mkdir -p TMP

${bed.name.endsWith(".gz")?"gunzip -c ":"cat"} "${bed}" |\\
	cut -f1,2,3,${column} |\\
	tail -n +2 > TMP/jeter.bed



if ${scatter_N_bed?true:false}
then
	sort -T TMP -t '\t' -k1,1 -k2,2n -S ${task.memory.kilo} TMP/jeter.bed > TMP/jeter.a

	sort -T TMP -t '\t' -k1,1 -k2,2n -S ${task.memory.kilo} "${scatter_N_bed}" > TMP/jeter.b

	bedtools subtract -sorted  -a TMP/jeter.a -b TMP/jeter.b |\\
		awk -F '\t' '(int(\$2) < int(\$3))' |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n -S ${task.memory.kilo} > TMP/jeter.bed
fi

## DELETION

awk -F '\t' '( \$4 <= (0.5 + ${treshold} ) )'  TMP/jeter.bed  |\\
	sort -T TMP -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n  |\\
	bedtools merge -d '${merge_d}' -o mean -c 4|\\
	awk -F '\t' '{
		L=int(\$3)-int(\$2)+1;
		if(L< ${min_size}) next;
		T=(\$4 <= ${treshold} ?"HOM_DEL":"HET_DEL");
		printf("%s\t%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,\$4,T);
		}' > TMP/del.bed

## DUP

awk -F '\t' '( \$4 >= (1.5 - ${treshold} ) )'  TMP/jeter.bed  |\\
	sort -T TMP -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n  |\\
	bedtools merge -d '${merge_d}' -o mean -c 4|\\
	awk -F '\t' '{
		L=int(\$3)-int(\$2)+1;
		if(L< ${min_size}) next;
		T=(\$4 >= 2.0 - ${treshold} ?"HOM_DUP":"HET_DUP");
		printf("%s\t%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,\$4,T);
		}' > TMP/dup.bed


cat TMP/del.bed TMP/dup.bed |\\
	sort -T TMP -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n > "${prefix}.bed"


cat << EOF > versions.yml
"${task.process}"
      awk: todo
EOF
"""
stub:
	def prefix = task.ext.prefix?:"${meta.id}.extract"
"""
touch versions.yml "${prefix}.bed"
"""
}

process BUILD_SVG {
tag "${meta.id}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
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

circle.MID {
        opacity: 0.8;
        stroke:gray;
        stroke-width:0.2;
		stroke-dasharray: 2, 2;
        fill:none;
        }
circle.DEL {
        opacity: 0.8;
        stroke:gray;
        stroke-width:0.2;
		stroke-dasharray: 2, 2;
        fill:none;
        }
		
circle.DUP {
        opacity: 0.8;
        stroke:gray;
        stroke-width:0.2;
		stroke-dasharray: 2, 2;
        fill:none;
        }
path.HET_DEL {
		opacity: 0.8;
        fill:red;
		stroke:gray;
        stroke-width:0.2;
        }
path.HOM_DEL {
		opacity: 0.8;
        fill:red;
		stroke:gray;
        stroke-width:0.2;
        }
path.HET_DUP {
	opacity: 0.8;
	fill:cyan;
	stroke:gray;
	stroke-width:0.2;
	}
path.HOM_DUP {
	opacity: 0.8;
	fill:blue;
	stroke:gray;
	stroke-width:0.2;
	}

line.cross {
        stroke:gray;
        stroke-width:1;
        }

EOF
cat "${moduleDir}/../../../src/css/fai.css" >> jeter.css

echo '<defs xmlns="${xmlns}">' > jeter.defs.xml
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

"""
stub:
    def prefix = task.ext.prefix?:"${meta.id}.circular"
"""
touch versions.yml  ${prefix}.svg
"""
}
