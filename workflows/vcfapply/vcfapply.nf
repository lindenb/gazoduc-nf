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

def gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults().putReference()

gazoduc.build("vcf", "NO_FILE").
	desc(gazoduc.Gazoduc.DESC_VCF_OR_VCF_LIST).
	existingFile().
	required().
	put()

gazoduc.build("bed","NO_FILE").
	desc("limit to that bed. Optional. If not defined, run on all chromosome").
	put()

gazoduc.build("script","NO_FILE").
	desc("user bash script. Takes a VCF stream on stdin as input and MUST produce an indexed VCF file 'selection.bcf'. A directory TMP can be used to handle temporary files").
	required().
	existingFile().
	put()

gazoduc.build("bed_cluster_method"," --size '1mb' ").
	desc("If a bed was provided, how do we split the BED using jvarkit/bedcluster").
	notEmpty().
	put()


include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {runOnComplete;moduleLoad;getVersionCmd} from '../../modules/utils/functions.nf'
include {SIMPLE_PUBLISH_01} from '../../modules/utils/publish.simple.01.nf'
include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {BCFTOOLS_CONCAT_01} from '../../subworkflows/bcftools/bcftools.concat.01.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'

if( params.help ) {
    gazoduc.usage().
	name("Vcf Apply").
	desc("Apply a script to a VCF file. Output is a VCF file").
	print();
    exit 0
} else {
   gazoduc.validate();
}


workflow {
	ch1 = APPLY_VCF_01(params, params.reference, file(params.script), file(params.vcf),file(params.bed))
	html = VERSION_TO_HTML(params,ch1.version)


	pub_ch = Channel.empty().mix(ch1.vcf).mix(ch1.version).mix(html.html)
	SIMPLE_PUBLISH_01(params, pub_ch.collect())
	}

workflow APPLY_VCF_01 {
     take:
        meta /* meta */
        reference /* indexed fasta reference */
	userScript
        vcfs /* file containing the path to VCF files . One per line */
	bed

     main:

        version_ch = Channel.empty()

	scr_ch = SCRIPT_VERSION(meta,userScript)
	version_ch = version_ch.mix(scr_ch.version)


	vcf2bed_ch = VCF_TO_BED(meta ,vcfs)
	version_ch = version_ch.mix( vcf2bed_ch.version)


	cluster_ch = CLUSTER(meta, reference, vcf2bed_ch. , bed)
	version_ch = version_ch.mix( cluster_ch.version)

	perctg_ch = PER_BED(meta, reference, userScript, cluster_ch.output.splitCsv(header:false,sep:"\t") )
	version_ch = version_ch.mix(perctg_ch.version)


	x3_ch = COLLECT_TO_FILE_01([:], perctg_ch.vcf.collect())
	version_ch = version_ch.mix(x3_ch.version)


	concat_ch = BCFTOOLS_CONCAT_01([:], x3_ch.output)
	version_ch = version_ch.mix(concat_ch.version)
		
		
	version_ch = MERGE_VERSION(meta, "vcffilterjdk", "vcffilterjdk",version_ch.collect())

	emit:
		vcf = concat_ch.vcf
		index = concat_ch.index
		version= version_ch
	}


runOnComplete(workflow)

process SCRIPT_VERSION {
executor "local"
input:
	val(meta)
	path(userScript)
output:
	path("version.xml"),emit:version
script:
"""
## https://stackoverflow.com/questions/12873682
echo -n "<properties><entry key='userScript'>" > version.xml

sed 's/&/\\&amp;/g; s/</\\&lt;/g; s/>/\\&gt;/g; s/"/\\&quot;/g; s/'"'"'/\\&#39;/g'  "${userScript}" >> version.xml

echo "</entry></properties>" >> version.xml
"""
}

process CLUSTER {
tag "${vcfbed}"
memory "3g"
input:
	val(meta)
	val(reference)
	path(vcfbed)
	path(bed)
output:
	path("beds.list"),emit:output
	path("version.xml"),emit:version
script:
	def method = meta.bed_cluster_method
"""
hostname 1>&2
${moduleLoad("bedtools bcftools jvarkit")}
set -o pipefail

mkdir -p BEDS TMP

cut -f1,2,3 "${vcfbed}" | sort -T TMP -t '\t' -k1,1 -k2,2n | bedtools merge > TMP/jeter.a.bed

if ${!bed.name.equals("NO_FILE")}

	${bed.name.endsWith(".gz")?"gunzip -c":"cat"} "${bed}" | cut  -f1,2,3 | sort -T TMP -t '\t' -k1,1 -k2,2n | bedtools merge > TMP/jeter.b.bed

	bedtools intersect -a TMP/jeter.a.bed -b TMP/jeter.b.bed |\
	java  -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar \${JVARKIT_DIST}/bedcluster.jar \
		-R "${reference}" \
		${method} \
		-o BEDS
else

	split --additional-suffix=.bed --lines=1 TMP/jeter.a.bed BEDS/cluster.
end

find ${PWD}/BEDS -type f -name "*.bed" > beds.list
test -s beds.list

##############################################################################################
cat <<- EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">Intersection with BED.</entry>
	<entry key="version">${getVersionCmd("bcftools bedtools")}</entry>
</properties>
EOF
"""
}


process PER_BED {
tag "${file(bed).name} ${file(vcf).name}"
memory "3g"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
	path(userScript)
	tuple val(bed),val(vcf)
output:
	path("selection.bcf"),optional:true,emit:vcf
	path("selection.bcf.csi"),optional:true,emit:index
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
mkdir -p BEDS TMP
${moduleLoad("bedtools bcftools")}
set -o pipefail

sort -T TMP -t '\t' -k1,1 -k2,2n "${vcfbed}" > TMP/jeter.a.bed
cut -f 1,2,3 "${bed}" | sort -T TMP -t '\t' -k1,1 -k2,2n | bedtools merge > TMP/jeter.b.bed

bedtools intersect -a TMP/jeter.a.bed -b TMP/jeter.b.bed > TMP/jeter.c.bed
cut -f 4 TMP/jeter.c.bed | uniq | sort -T TMP | uniq >  TMP/vcfs.list
cut -f1,2,3 TMP/jeter.c.bed | sort -T TMP -t '\t' -k1,1 -k2,2n | bedtools merge > TMP/jeter.d.bed

if test  -s TMP/vcfs.list ; then

	bcftools concat --regions-file "TMP/jeter.d.bed" -a --file-list TMP/vcfs.list -O u | bash ${userScript.toRealPath()}  "${reference}"
fi


##############################################################################################
cat <<- EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">Apply user script</entry>
	<entry key="version">${getVersionCmd("bcftools bedtools")}</entry>
</properties>
EOF
"""
}
