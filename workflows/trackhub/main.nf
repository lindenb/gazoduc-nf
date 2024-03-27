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


include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {escapeXml;runOnComplete;moduleLoad;dumpParams} from '../../modules/utils/functions.nf'


if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}


workflow {
	ch1 = TRACK_HUB(params.genomeId, Channel.fromPath(params.vcfs))
	}

runOnComplete(workflow)

workflow TRACK_HUB {
	take:
		genomeId
		vcfs
	main:
		version_ch = Channel.empty()

		concat0_ch= CONCAT_VCFS1(vcfs.splitText().map{it.trim()}.collate(200))
		concat_ch= CONCAT_VCFS2(concat0_ch.output.collect())

		c1_ch = COMPILE_JAVA()
                version_ch = version_ch.mix(c1_ch.version)
		
		ci_ch = DOWNLOAD_CHROM_INFO(genomeId)
		version_ch = version_ch.mix(ci_ch.version)

		scan_ch = SCAN_VCF( c1_ch.output, concat_ch.output )
		version_ch = version_ch.mix(scan_ch.version)

		each_type = scan_ch.types.splitText().map{it.trim()}

		all_outputs = Channel.empty()
		bed_ch = MERGE_BED(ci_ch.output.combine(vcfs.combine(scan_ch.bed.combine(each_type))) )
		version_ch = version_ch.mix(bed_ch.version)
		all_outputs = all_outputs.mix(bed_ch.output)
	
		if ( params.with_bnd as boolean ) {
			inter_ch = MERGE_BED_PE(ci_ch.output , vcfs, scan_ch.bedpe)
			version_ch = version_ch.mix(inter_ch.version)
			all_outputs = all_outputs.mix(inter_ch.output)
			}

		track_ch = MAKE_TRACK(genomeId, all_outputs.collect() )
		version_ch = version_ch.mix(track_ch.version)

		version_ch = MERGE_VERSION("Trackhub", version_ch.collect())
	emit:
		version = version_ch
		tar = track_ch.tar
	}


process CONCAT_VCFS1 {
tag "N=${L.size()}"
afterScript "rm -rf TMP"
memory '5G'
input:
	val(L)
output:
	path("concat.vcf.gz"),emit:output
script:
"""
hostname 1>&2
${moduleLoad("jvarkit bcftools")}
mkdir -p TMP

cat << EOF > TMP/jeter.list
${L.join("\n")}
EOF

java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfconcat --merge	< TMP/jeter.list |\\
	bcftools view -O z -o TMP/jeter.vcf.gz

mv -v TMP/jeter.vcf.gz concat.vcf.gz
"""
}


process CONCAT_VCFS2 {
tag "N=${L.size()}"
afterScript "rm -rf TMP"
memory '5G'
input:
	val(L)
output:
	path("concat.bcf"),emit:output
script:
"""
hostname 1>&2
${moduleLoad("jvarkit bcftools")}
mkdir -p TMP


cat << EOF > TMP/jeter.list
${L.join("\n")}
EOF

java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfconcat --merge	< TMP/jeter.list |\\
	bcftools view -O u -o TMP/jeter.bcf

bcftools sort -T TMP/sort --max-mem ${task.memory.giga}G -O b -o TMP/jeter2.bcf TMP/jeter.bcf

mv -v TMP/jeter2.bcf concat.bcf

"""
}


process COMPILE_JAVA {
	afterScript "rm -rf TMP"
	output:
		path("minikit.jar"),emit:output
		path("version.xml"),emit:version
	script:
	"""
	hostname 1>&2
	${moduleLoad("jvarkit")}
	
	mkdir -p TMP

	cat "${moduleDir}/Minikit.java" |\\
		sed 's/__MIN_SV_LEN__/${params.min_sv_len}/;s/__MAX_SV_LEN__/${params.max_sv_len}/;s/__WITH_BND__/${params.with_bnd as boolean}/' >  TMP/Minikit.java

	javac -d TMP -cp \${JVARKIT_DIST}/jvarkit.jar -sourcepath TMP  TMP/Minikit.java
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


process SCAN_VCF {
	tag "${vcf.name}"
	afterScript "rm -rf TMP"
	input:
		path(minikit)
		path(vcf)
	output:
		path("all.bed"),emit:bed
		path("all.bedpe"),emit:bedpe
		path("types.txt"),emit:types
		path("version.xml"),emit:version
	when:
		true
	script:
	"""
	hostname 1>&2
	${moduleLoad("bcftools jvarkit")}
	mkdir -p TMP
	
	bcftools view "${vcf}" |\
		java -Djava.io.tmpdir=TMP -cp \${JVARKIT_DIST}/jvarkit.jar:${minikit} Minikit "${vcf}" TMP/tmp.bed TMP/tmp.bedpe

	cut -f 10 TMP/tmp.bed | uniq | LC_ALL=C sort -T TMP | LC_ALL=C uniq > types.txt

	mv -v TMP/tmp.bed ./all.bed
	mv -v TMP/tmp.bedpe ./all.bedpe

	###############################################################################
	cat << EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">extract variants</entry>
	</properties>
	EOF
	"""
	}




process DOWNLOAD_CHROM_INFO {
executor "local"
input:
	val(genomeId)
output:
	path("chromInfo.txt"),emit:output
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def url = "https://hgdownload.cse.ucsc.edu/goldenpath/${genome.ucsc_name}/database/chromInfo.txt.gz"
"""
wget -O "chromInfo.txt.gz" "${url}"
gunzip chromInfo.txt.gz

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">download chrominfo</entry>
	<entry key="url">${url}</entry>
</properties>
EOF
"""
}


process MERGE_BED {
	tag "${type} ${bed.name}"
	afterScript "rm -rf TMP"
	input:
		tuple path(chrominfo),path(vcfs),path(bed),val(type)
	output:
		path("output.tsv"),emit:output
		path("version.xml"),emit:version
	script:
		def filename = file("${params.prefix?:""}sv.${type}.bb")
		def now =  new java.text.SimpleDateFormat("yyyy-MM-dd").format(new java.util.Date());

	"""
	hostname 1>&2
	${moduleLoad("ucsc")}
	set -o pipefail
	mkdir -p TMP

cat << EOF > TMP/schema.as
table sv
"Structural Variant"
    (
    string chrom;      "Chromosome"
    uint chromStart;   "Start position"
    uint chromEnd;     "End position"
    string name;       "Sample Name"
    uint score;        "Score from 0-1000."
    char[1] strand;    "+ or -"
    uint thickStart;   "VCF/POS"
    uint thickEnd;   "VCF/END"
    uint reserved;        "rgb color of item"
    string type; "Sv type=${type}"
    int svLen; "Sv Len"
    string filters; "FILTERs"
    string gtType; "Genotype Type"
    string vcfName; "VCF source"
    )
EOF


	awk -F '\t' '(\$10=="${type}")' '${bed}' |\\
	LC_ALL=C sort -T TMP -k1,1 -k2,2n > TMP/all.bed

	bedToBigBed -as=TMP/schema.as -type=bed9+5 TMP/all.bed "${chrominfo}" "${filename.name}"



#################################################

cat << EOF > trackDb.txt
track ${filename.getBaseName()}
bigDataUrl ${filename.name}
shortLabel ${filename.getBaseName()}
longLabel ${filename.getBaseName()} ${type} (${params.min_sv_len} < SVLEN < ${params.max_sv_len})
type bigBed 9
itemRgb "On"

EOF

#################################################

cat << EOF > TMP/jeter.html
<html>
<body>
<h1>${filename.name}</h1>
<p>${type} (${params.min_sv_len} < SVLEN < ${params.max_sv_len}</p>
<p>${escapeXml(params.description)}</p>
<h3>Schema</h3><pre>
EOF

cat TMP/schema.as >> TMP/jeter.html

cat << EOF >> TMP/jeter.html
</pre><h3>Sources</h3>
<pre>
EOF

cat '${vcfs}' >> TMP/jeter.html

cat << EOF >> TMP/jeter.html
</pre>
<hr/>
<div>Author: Pierre Lindenbaum. Generated on ${now}</div>
</body>
</html>
EOF

mv -v TMP/jeter.html "${filename.getBaseName()}.html"



	echo "\${PWD}/trackDb.txt\t\${PWD}/${filename.getBaseName()}.html\t\${PWD}/${filename.name}" > output.tsv


	###############################################################################
	cat << EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="type">${type}</entry>
		<entry key="description">create SV.as</entry>
	</properties>
	EOF
	"""
	}




process MERGE_BED_PE {
	tag "${bed.name}"
	afterScript "rm -rf TMP"
	input:
		path(chrominfo)
		path(vcfs)
		path(bed)
	output:
		path("output.tsv"),emit:output
		path("version.xml"),emit:version
	script:
		def filename = file("${params.prefix?:""}interact.bb")
		def now =  new java.text.SimpleDateFormat("yyyy-MM-dd").format(new java.util.Date());
	"""
	hostname 1>&2
	${moduleLoad("ucsc/0.0.0")}
	mkdir -p TMP


cat << EOF > TMP/schema.as
table interact
"BND"
    (
    string chrom;      "Chromosome (or contig, scaffold, etc.). For interchromosomal, use 2 records"
    uint chromStart;   "Start position of lower region. For interchromosomal, set to chromStart of this region"
    uint chromEnd;     "End position of upper region. For interchromosomal, set to chromEnd of this region"
    string name;       "Name of item, for display.  Usually 'sourceName/targetName' or empty"
    uint score;        "Score from 0-1000."
    double value;      "Strength of interaction or other data value. Typically basis for score"
    string exp;        "Experiment name (metadata for filtering). Use . if not applicable"
    string color;      "Item color.  Specified as r,g,b or hexadecimal #RRGGBB or html color name, as in //www.w3.org/TR/css3-color/#html4."
    string sourceChrom;  "Chromosome of source region (directional) or lower region. For non-directional interchromosomal, chrom of this region."
    uint sourceStart;  "Start position source/lower/this region"
    uint sourceEnd;    "End position in chromosome of source/lower/this region"
    string sourceName;  "Identifier of source/lower/this region"
    string sourceStrand; "Orientation of source/lower/this region: + or -.  Use . if not applicable"
    string targetChrom; "Chromosome of target region (directional) or upper region. For non-directional interchromosomal, chrom of other region"
    uint targetStart;  "Start position in chromosome of target/upper/this region"
    uint targetEnd;    "End position in chromosome of target/upper/this region"
    string targetName; "Identifier of target/upper/this region"
    string targetStrand; "Orientation of target/upper/this region: + or -.  Use . if not applicable"
    string filters; "FILTER column"
    string gtType; "Genotype Type"
    string vcfName; "VCF source"
    )
EOF


	LC_ALL=C sort -T TMP -k1,1 -k2,2n "${bed}" | uniq > TMP/all.bedpe

	bedToBigBed -as=TMP/schema.as -type=bed5+16 TMP/all.bedpe "${chrominfo}"  "${params.prefix?:""}interact.bb"



cat << EOF > trackDb.txt

track ${filename.getBaseName()}
bigDataUrl ${filename.name}
shortLabel ${filename.getBaseName()}
longLabel ${filename.getBaseName()} BND
type bigInteract
itemRgb "On"

EOF


#################################################

cat << EOF > TMP/jeter.html
<html>
<h1>${filename.name}</h1>
<p>BND</p>
<p>${escapeXml(params.description)}</p>
<h3>Schema</h3><pre>
EOF

cat "TMP/schema.as" >> TMP/jeter.html

cat << EOF >> TMP/jeter.html
</pre><h3>Sources</h3>
<pre>
EOF

cat '${vcfs}' >> TMP/jeter.html

cat << EOF >> TMP/jeter.html
</pre>
<hr/>
<div>Author: Pierre Lindenbaum. Generated on ${now}</div>
</html>
EOF


mv -v "TMP/jeter.html" "${filename.getBaseName()}.html"


	echo "\${PWD}/trackDb.txt\t\${PWD}/${filename.getBaseName()}.html\t\${PWD}/${filename.name}" > output.tsv


	###############################################################################
	cat << EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
	</properties>
	EOF
	"""
	}

process MAKE_TRACK {
	tag "N=${L.size()}"
	//afterScript "rm -rf TMP"
	input:
		val(genomeId)
		val(L)
	output:
		path("${params.prefix?:""}hub.tar.gz"),emit:tar
		path("version.xml"),emit:version
	script:
		def build = params.genomes[genomeId].ucsc_name
		def prefix= params.prefix?:""
		def now =  new java.text.SimpleDateFormat("yyyy-MM-dd").format(new java.util.Date());
	"""
	hostname 1>&2
	${moduleLoad("ucsc/0.0.0")}
	mkdir -p "TMP/${build}"
	set -x

	cat ${L.join(" ")} > TMP/jeter.list

#################################################

cat << EOF > TMP/genomes.txt
genome ${build}
trackDb ${build}/trackDb.txt
EOF

#################################################

cat << EOF > TMP/hub.html
<html>
<head><title>${prefix}hub</title></head>
<body>
<h1>${prefix}hub</h1>
<p>${escapeXml(params.description)}</p>
<hr/>
<div>Author: Pierre Lindenbaum. Generated on ${now}</div>
</body>
</html>
EOF

#################################################

cat << EOF > TMP/hub.txt
hub ${prefix}hub
shortLabel ${prefix}hub
longLabel ${prefix}hub
genomesFile genomes.txt
email plindenbaum@yahoo.fr
descriptionUrl hub.html 
EOF

touch TMP/${build}/trackDb.txt

cut -f 1 TMP/jeter.list | xargs -L 50 cat  >> TMP/${build}/trackDb.txt

cut -f 2 TMP/jeter.list	| xargs -L1 -I '{}' cp -v '{}' TMP/${build}/
cut -f 3 TMP/jeter.list	| xargs -L1 -I '{}' cp -v '{}' TMP/${build}/

#################################################

find TMP -type f 1>&2

hubCheck "\${PWD}/TMP/hub.txt"

rm TMP/jeter.list
mv -v TMP "${prefix}hub"

tar cvfz "${prefix}hub.tar.gz" "${prefix}hub"

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
</properties>
EOF
"""
}
