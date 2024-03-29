/*

Copyright (c) 2023 Pierre Lindenbaum

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
include {getVersionCmd;moduleLoad;runOnComplete;dumpParams} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {SAMTOOLS_SAMPLES} from '../../subworkflows/samtools/samtools.samples.03.nf'


if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}

workflow {
	PLOT_WGS_COVERAGE(params, file(params.bams))
	}

runOnComplete(workflow)

workflow PLOT_WGS_COVERAGE {
	take:
		meta
		bams
	main:
		version_ch = Channel.empty()
		to_zip = Channel.empty()

		snbam_ch = SAMTOOLS_SAMPLES(["with_header":true,"allow_multiple_references":true,"allow_duplicate_samples":true],  bams)
		version_ch = version_ch.mix(snbam_ch.version)

		plot_ch = PLOT_WGS([:],snbam_ch.rows)
		version_ch = version_ch.mix(plot_ch.version)
		to_zip = to_zip.mix(plot_ch.svg)

		version_ch = MERGE_VERSION([:], "wgsplot", "WGS plot", version_ch.collect())
		to_zip = to_zip.mix(version_ch)

		html =  VERSION_TO_HTML(version_ch)
		to_zip = to_zip.mix(html.html)

		zip_ch = ZIPIT([:],to_zip.collect())
	emit:
		version = version_ch
		zip = zip_ch.zip
	}



process PLOT_WGS {
tag "${row.new_sample}"
memory '5g'
afterScript "rm -f jeter.svg"
input:
        val(meta)
	val(row)
output:
	path("${params.prefix?:""}${row.new_sample}.svg"),emit:svg
	path("version.xml"),emit:version
script:
	def maxcov = params.maxcov?:80
	def mapq = params.mapq?:80
	def dimension = params.dimension?:"3000x400"
"""
	hostname 1>&2
	${moduleLoad("jvarkit")}
	set -x

	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar \${JVARKIT_DIST}/wgscoverageplotter.jar \
		--mapq ${mapq} \
		${params.include_contig_regex.isEmpty()?"":"--include-contig-regex \"${params.include_contig_regex}\""} \
		--dimension "${dimension}" \
		-C "${maxcov}" -R "${row.reference}" "${row.bam}" > jeter.svg

	mv jeter.svg "${params.prefix?:""}${row.new_sample}.svg"

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">wgs coverage plotter</entry>
	<entry key="sample">${row.new_sample}</entry>
	<entry key="maxcov">${maxcov}</entry>
	<entry key="mapq">${mapq}</entry>
	<entry key="bam">${row.bam}</entry>
	<entry key="version">${getVersionCmd("jvarkit/wgscoverageplotter")}</entry>
</properties>
EOF
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
