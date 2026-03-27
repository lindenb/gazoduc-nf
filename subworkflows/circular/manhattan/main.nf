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
include { FAI_TO_SVG            } from '../../../subworkflows/circular/fai2svg'
include { CYTOBAND_TO_SVG       } from '../../../subworkflows/circular/cytoband2svg'

/** 
   plot each bed as a single circular manhattan plot */
workflow CIRCULAR_MANHATTAN {
take:
    metadata
    fasta
    fai
    dict
    gtf

    beds //meta,bed
main:
    versions = Channel.empty()
    multiqc  = Channel.empty()
    
    FAI_TO_SVG(metadata,fai,dict)
    versions = versions.mix(FAI_TO_SVG.out.versions)


    CYTOBAND_TO_SVG(
        metadata,
        FAI_TO_SVG.out.dict,
        Channel.empty() //empty cytoband, will ve downloaded
        )
    versions = versions.mix(CYTOBAND_TO_SVG.out.versions)

    BED_TO_XML(
        FAI_TO_SVG.out.dict.first(),
        beds
        )
    versions = versions.mix(BED_TO_XML.out.versions)

    XSLTPROC(
            [[id:"manhattan2svg"],"${moduleDir}/../../../src/xsl/manhattan2svg.xsl"],
            BED_TO_XML.out.xml
            )
    versions = versions.mix(XSLTPROC.out.versions)

    BUILD_SVG(
        FAI_TO_SVG.out.svg.first(),
        CYTOBAND_TO_SVG.out.svg.first(),
        XSLTPROC.out.xml
        )
    versions = versions.mix(BUILD_SVG.out.versions)

emit:
    versions
    multiqc
    svg = BUILD_SVG.out.svg
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

text.title {
        font-size:12px;
        text-anchor: middle;
        }



circle { 
        fill:gray;  
        stroke:darkgray;
        opacity:0.4;
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
