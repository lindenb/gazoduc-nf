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

/** workflow is just for quick dev/testing */
include {DOWNLOAD_CYTOBAND    } from '../../modules/ucsc/download.cytobands'
include {FAI2BED              } from '../../modules/samtools/fai2bed'
include {DOWNLOAD_GTF_OR_GFF3 as DOWNLOAD_GTF} from '../../modules/gtf/download'
include {GTF_TO_BED           } from '../../modules/jvarkit/gtf2bed'


jvarkit = "java -jar \${HOME}/jvarkit.jar"

workflow {
    fasta=[[id:"GRCH38",ucsc_name:"hg38"],params.fasta]
    fai = [ fasta[0], params.fasta+".fai"]
    dict= [ fasta[0], params.dict]
    
    DOWNLOAD_CYTOBAND(dict)
    FAI2BED(fai)
    DOWNLOAD_GTF(dict)
    GTF_TO_BED(dict,DOWNLOAD_GTF.out.gtf.map{meta,gtf,tbi->[meta,gtf]})

    FAI_TO_SVG(dict, FAI2BED.out.bed)
    CYTOBANDS_TO_SVG(FAI_TO_SVG.out.dict, DOWNLOAD_CYTOBAND.out.bed)
    GTF_TO_SVG(FAI_TO_SVG.out.dict, GTF_TO_BED.out.bed)

    vcfs_ch = Channel.fromPath(params.vcfs).splitText()
	.map{it.trim()}
	.collect(sort:true)
	.flatMap{v->
		def L=[];
		for(int i=0;i<v.size();i++) {
			L.add ( [ [id:v[i].md5(),index:i,count:v.size()],file(v[i])] );
			}
		return L;
		}
	
    
    VCF_TO_SVG(FAI_TO_SVG.out.dict, vcfs_ch)


    MAKE_SVG(
        FAI_TO_SVG.out.svg,
        CYTOBANDS_TO_SVG.out.svg,
        GTF_TO_SVG.out.svg,
	VCF_TO_SVG.out.svg.map{meta,svg->svg}.collect().map{[[id:"cnvs"],it.sort()]},
	VCF_TO_SVG.out.css.first()
        )
}

process FAI_TO_SVG {
    executor "local"
    //cache false
    input:
        tuple val(meta1),path(dict)
        tuple val(meta ),path(bed)
    output:
        tuple val(meta),path("*.fai.svg"),path("*.css"),path("*.defs.svg"), emit:svg
        tuple val(meta),path("*.dict"),emit:dict
        path("versions.yml"),emit:versions
    script:
        def stylesheet= "${moduleDir}/../../src/xsl/fai2svg.xsl"
        def css = "${moduleDir}/../../src/css/fai.css"

    """
    ${jvarkit} bed2xml  --min-contig-length 50000000 --dict-out ${meta.id}.out.dict -R ${dict} "${bed}" -o jeter.xml
    xmllint  --noout jeter.xml
    xmllint  --noout ${stylesheet}
    xsltproc --stringparam radius 1000  -o "${meta.id}.fai.svg"  ${stylesheet} jeter.xml
    
    touch ${meta.id}.defs.svg
    cp "${css}" "${meta.id}.css"
    touch versions.yml
    """
    }




process CYTOBANDS_TO_SVG {
     executor "local"
    //cache false
    input:
        tuple val(meta1),path(dict)
        tuple val(meta ),path(bed)
    output:
        tuple val(meta),path("cytoband.svg"),path("cytoband.css"),path("cytoband.defs.svg"), emit:svg
        path("versions.yml"),emit:versions
    script:
        def stylesheet= "${moduleDir}/../../src/xsl/cytoband2svg.xsl"

    """
    ${jvarkit} bed2xml --columns "name,color" -R ${dict} "${bed}" -o jeter.xml
    xmllint  --noout jeter.xml
    xmllint  --noout ${stylesheet}
    xsltproc \\
        --stringparam radius_R1 1000 \\
        --stringparam radius_R2 950 \\
         -o "cytoband.svg"  ${stylesheet} jeter.xml
    
    # cache1
    touch cytoband.defs.svg
    touch cytoband.css
    touch versions.yml
    """
    }


process GTF_TO_SVG {
    executor "local"
    //cache false
    input:
        tuple val(meta1),path(dict)
        tuple val(meta ),path(bed)
    output:
        tuple val(meta),path("gtf.svg"),path("gtf.css"),path("gtf.defs.svg"), emit:svg
        path("versions.yml"),emit:versions
    script:
        def stylesheet= "${moduleDir}/../../src/xsl/gtf2svg.xsl"

    """
    ${jvarkit} bed2xml  --columns gene_name,gene_biotype --distance 10000 -R ${dict} "${bed}" -o jeter.xml
   xmllint  --noout jeter.xml
   xmllint  --noout ${stylesheet}
   xsltproc \\
        --stringparam radius_R1 900 \\
        --stringparam radius_R2 850 \\
         -o "gtf.svg"  ${stylesheet} jeter.xml


    # cache1

cat << EOF > gtf.css
path.gene {
    stroke:none;
    fill: lightgray;
    }

path.protein_coding {
    fill: green;
    }

path.lncRNA {
    fill: blue;
    }

EOF

    touch gtf.defs.svg
    touch versions.yml
    """
    }



process VCF_TO_SVG {
	tag "${meta.id}"
    executor "local"
    //cache false
    conda "${moduleDir}/../../conda/bioinfo.01.yml"
    input:
        tuple val(meta1),path(dict)
        tuple val(meta ),path(vcf)
    output:
        tuple val(meta),path("*.svg"),emit:svg
	tuple val(meta),path("vcf.css"),emit:css
        path("versions.yml"),emit:versions
    script:
        def stylesheet= "${moduleDir}/../../src/xsl/dragen.sv.xsl"

	def index = (meta.index as int)
	def count = (meta.count as int)
	def r0 = 850.0
	def r1 = r0 * 0.1
	def dr = (r0-r1)/count
	def R0 = r0 - (dr * index)
	def R1 = R0 - Math.min(50.0,dr*0.95)
    def fill = (index%2==0?"gainsboro":"lavenderBlush")
    """
    NAME=\$(bcftools query -l "${vcf}")
    bcftools view -O u --apply-filters '.,PASS' -m2 '${vcf}' |\\
    bcftools query -f '%CHROM\t%POS0\t%END\t%ALT\n' |\\
	tr -d '<>' |\\
    	${jvarkit} bed2xml  --columns svtype --distance 10000 -R ${dict}  -o jeter.xml
   xmllint  --noout jeter.xml
   xmllint  --noout ${stylesheet}
   xsltproc \\
	--stringparam title "\${NAME}" \\
        --stringparam radius_R1 ${R0} \\
        --stringparam radius_R2 ${R1} \\
        --stringparam style "opacity:0.7;stroke-dasharray:0 4 0;stroke:silver;fill:${fill};" \\
         -o "${meta.id}.svg"  ${stylesheet} jeter.xml


    # cache9

cat << EOF > vcf.css


path.cnv {
    stroke: blue;
    stroke-width: 2px;
    fill: lightgray;
    }

path.DUP {
    fill: green;
    }

path.DEL {
    fill: blue;
    }

EOF

    touch versions.yml
    """
    }



process MAKE_SVG {
    executor "local"
    cache false
    input:
        tuple val(meta),path(fai_svg),path(fai_css),path(fai_defs)
        tuple val(meta2),path(cytoband_svg),path(cytoband_css),path(cytoband_defs)
        tuple val(meta3),path(gtf_svg),path(gtf_css),path(gtf_defs)
	tuple val(meta ),path("SVGS/*")
	tuple val(meta4),path(xcss)
    output:
        tuple val(meta),path("*.svg"),emit:svg
        path("versions.yml"),emit:versions
    script:
        def xmlns = "http://www.w3.org/2000/svg"
        def margin = (task.ext.margin?:100) as int
        def margin_top = (task.ext.margin_top?:margin) as int
        def margin_left = (task.ext.margin_left?:margin) as int
        def margin_bottom = (task.ext.margin_bottom?:margin) as int
        def margin_right = (task.ext.margin_right?:margin) as int
        def radius = (task.ext.radius?:1000) as int
        def diameter = (radius*2)
        def image_width = diameter + margin_left + margin_right
        def image_height = diameter + margin_top + margin_bottom
        def title = task.ext.title?:"${meta.title?:meta.id}"
        def prefix = task.ext.prefix?:"${meta.id}"
        
    """
mkdir -p TMP


cat SVGS/*.svg  | awk 'BEGIN {print("<g>");} {print;} END{print("</g>");}' > extra.svg

cat ${fai_defs} ${cytoband_defs} ${gtf_defs} | awk 'BEGIN {print("<defs>");} {print;} END{print("</defs>");}' > defs.svg

cat << EOF > style.css
rect.bckg {
    fill:none;
    stroke:gray;
    }
EOF

cat ${fai_css} ${cytoband_css} ${gtf_css} ${xcss} >> style.css

cat << __EOF__ > jeter.svg
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
  <xi:include parse="text" href="style.css"/>
  </style>
  <xi:include href="defs.svg"/>
  <g>
    <rect x="0" y="0" width="${image_width +1}" height="${image_width +1}" class="bckg"/>
    <g transform="translate(${margin_left},${margin_top})">
        <g transform="translate(${radius},${radius})">
            <xi:include href="${fai_svg}"/>
            <xi:include href="${cytoband_svg}"/>
            <xi:include href="${gtf_svg}"/>
            <xi:include href="extra.svg"/>
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


xmllint --format --xinclude jeter.svg > jeter2.svg
mv jeter2.svg jeter.svg


mv jeter.svg ${prefix}.svg
cat ${prefix}.svg  |gzip --best > \${HOME}/jeter.svg.gz
touch versions.yml
"""
}
