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
include {isGRCh37;isGRCh38} from '../../../modules/utils/k1.nf'


workflow REVEL {
	take:
		meta
		fasta
		fai
		dict
		vcfs /* meta,vcf,idx */
	main:
		versions = Channel.empty()
		
		if(isGRCh37(fai[1])) {
			column =  2
		}else if(isGRCh38(fai[1])) {
			column =  3
		} else {
			column =  -1
		}
		
		DOWNLOAD(fasta,fai,dict,column)
		versions = versions.mix(DOWNLOAD.out.versions)
		
		ANNOTATE(DOWNLOAD.out.output, vcfs)
		versions = versions.mix(ANNOTATE.out.versions)
	emit:
		vcf = ANNOTATE.out.vcf
		versions
}

process DOWNLOAD {
tag "${meta1.id?:fasta.name}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	val(colpos)
output:
	tuple val(meta1),path("*.bed.gz"),path("*.bed.gz.tbi"),path("*.header"),emit:output
	path("versions.yml"),emit:versions
script:
	def TAG = "REVEL"
	def WHATIZ="REVEL is an ensemble method for predicting the pathogenicity of missense variants in the human genome. https://doi.org/10.1016/j.ajhg.2016.08.016 ."
	def url = task.ext.url?:"https://zenodo.org/records/7072866/files/revel-v1.3_all_chromosomes.zip?download=1"
	//def colpos= isGRCh37(fai)?2: isGRCh38(fai)?3:-1
	if(colpos<1) {
		throw new IllegalArgumentException("Bad fai version ${task.process}");
		}
"""
hostname 1>&2
mkdir -p TMP
wget -O TMP/jeter.zip "${url}"


unzip -p TMP/jeter.zip  revel_with_transcript_ids |\\
	tail -n +2 |\\
	tr "," "\t" |\\
	cut -f 1,${colpos},4,5,8 |\\
	awk -F '\t' '\$2!="."' |\\
	uniq |\\
	jvarkit  -Djava.io.tmpdir=TMP  bedrenamechr -R "${fasta}" --column 1 --convert SKIP  |\\
		LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n |\\
		uniq |\\
		bgzip > TMP/${TAG}.bed.gz

tabix -s 1 -b 2 -e 2 -f TMP/${TAG}.bed.gz


mv TMP/${TAG}.bed.gz ./
mv TMP/${TAG}.bed.gz.tbi ./

echo '##INFO=<ID=${TAG},Number=1,Type=Float,Description="${WHATIZ} ${url}">' > ${TAG}.header


cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${url}"
END_VERSIONS
"""
}



process ANNOTATE {
tag "${meta.id?:vcf.name}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(tabix),path(tbi),path(header)
	tuple val(meta ),path(vcf),path(vcf_idx)
output:
	tuple val(meta),path("*.bcf"),path("*.bcf.csi"),emit:vcf
	path("versions.yml"),emit:versions
script:
	def TAG = "REVEL"
	def prefix = task.ext.prefix?:vcf.baseName+".revel"
"""
mkdir -p TMP OUTPUT

bcftools annotate \\
	--threads ${task.cpus} \\
	-a "${tabix}" \\
	-h "${header}" \\
	-c "CHROM,POS,REF,ALT,${TAG}" \\
	-O b \\
	--merge-logic '${TAG}:max' \\
	-o TMP/${prefix}.bcf \\
	'${vcf}'

bcftools index \\
	-f \\
	--threads ${task.cpus} \\
	TMP/${prefix}.bcf

mv TMP/${prefix}.bcf ./
mv TMP/${prefix}.bcf.csi ./


cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
}
