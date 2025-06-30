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



workflow ANNOTATE_BHFUCL {
	take:
		meta
		fasta
		fai
		dict
		gtf
		vcfs /* meta, vcf,vcf_index */
	main:
		source_ch = DOWNLOAD(fasta,fai,dict,gtf)
		annotate_ch = ANNOTATE(source_ch.bed, source_ch.tbi,vcfs)
	emit:
		output = annotate_ch.output
		doc = source_ch.doc
	}
		
process DOWNLOAD {
tag "${fasta.name}"
afterScript "rm -rf TMP"
label "process_quick"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta4),path(gtf),path(gtf_tbi)
output:
	tuple val(meta1),path("*.bed.gz"),emit:bed
	tuple val(meta1),path("*.bed.gz.tbi"),emit:tbi
	tuple val(meta1),path("*.md"),emit:doc
script:
    def TAG = task.ext.tag?:"BHFUCL"
    def URL = "http://ftp.ebi.ac.uk/pub/databases/GO/goa/bhf-ucl/gene_association.goa_bhf-ucl.gz"
    def WHATIZ = "Cardiovascular Gene Ontology Annotation Initiative ${URL}"
"""
hostname 1>&2
mkdir -p TMP
set -o pipefail
mkdir -p TMP

jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP gtf2bed  --columns "gtf.feature,gene_name" -R "${fasta}"  "${gtf}" |\\
	awk -F '\t' '\$4=="gene" && \$5!="." && \$5!=""' |\\
	cut -f1,2,3,5 |\\
        LC_ALL=C sort --buffer-size=${task.memory.mega}M -t '\t' -k4,4 -T TMP  |\\
	uniq > TMP/genes.bed

wget -O - "${URL}" |\\
	gunzip -c |\\
	awk -F '\t' '/#/{next} (\$4 ~ /^NOT/) {next;} {print \$3;}' |\\
	uniq | LC_ALL=C sort -T TMP | uniq > TMP/genes.txt

LC_ALL=C  join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,1.4' TMP/genes.bed TMP/genes.txt |\\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\\
	uniq |\\
	bgzip > TMP/${TAG}.bed.gz


tabix --force -p bed TMP/${TAG}.bed.gz

mv TMP/${TAG}.bed.gz ./
mv TMP/${TAG}.bed.gz.tbi ./


cat << __EOF__ > ${TAG}.md
${WHATIZ}
__EOF__
"""
}



process ANNOTATE {
tag "${vcf.name}"
afterScript "rm -rf TMP"
label "process_quick"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(tabix)
	tuple val(meta2),path(tbi)
	tuple val(meta),path(vcf),path(vcf_idx)
output:
    tuple val(meta),path("*.bcf"),path("*.csi"),emit:output
script:
    def TAG = task.ext.tag?:"BHFUCL"
    def distance = 1000;
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP OUTPUT

bcftools view '${vcf}' |\\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${HOME}/packages/jvarkit/dist/jvarkit.jar  \
		vcfnearest  \\
			--bed '${tabix}' \\
			--tag "${TAG}_NEAR" \\
			--distance '${distance}' \\
			-C ',,,GENE' |\\
	bcftools view  -O b -o TMP/${TAG}.${vcf.getBaseName()}.bcf
    
bcftools index \\
    --threads ${task.cpus} \\
    --force TMP/${TAG}.${vcf.getBaseName()}.bcf

mv TMP/*.bcf ./
mv TMP/*.bcf.csi ./
"""
}
