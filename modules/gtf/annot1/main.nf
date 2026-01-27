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
//include { parseBoolean } from '../../utils/functions.nf'
//include { isBlank      } from '../../utils/functions.nf'


/*
Basic annotation of a VCF with a GTF (in_cds, in intron, etc...)
*/
process GTF_ANNOTATION {
label "process_single"
tag "${meta.id}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(fai)//for slop
	tuple val(meta2),path(gtf)
	tuple val(meta3),path(genes_of_interest) // optional
	tuple val(meta ),path(vcf)
output:
	tuple val(meta),path("*.vcf.gz"),emit:vcf
	tuple val(meta),path("*.roi.bed"),optional:true,emit:roi_bed
    path("versions.yml"),emit:versions
script:
	def xstream = task.ext.xstream?:"1000"
    def prefix = task.ext.prefix?:"${meta.id}.annot"
    def args1 = task.ext.args1?:""
    def args2 = task.ext.args2?:""
	def jvm =  task.ext.jvm?:" -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData"

"""
mkdir -p TMP

set -x

${gtf.name.endsWith(".gz")?"gunzip -c":"cat"} "${gtf}" |\\
	awk -F '\t' '(\$3=="exon") {printf("%s\t%d\t%s\\n",\$1,int(\$4)-1,int(\$5))}' |\\
	sort -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n -T TMP |\\
	bedtools merge  > TMP/exon.bed

echo '##INFO=<ID=in_exon,Number=0,Type=Flag,Description="In exon of ${gtf.name}">' > TMP/exon.hdr

${gtf.name.endsWith(".gz")?"gunzip -c":"cat"} "${gtf}" |\\
	awk -F '\t' '(\$3=="CDS") {printf("%s\t%d\t%s\\n",\$1,int(\$4)-1,int(\$5))}' |\\
	sort -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n -T TMP |\\
	bedtools merge > TMP/cds.bed

echo '##INFO=<ID=in_cds,Number=0,Type=Flag,Description="In CDS of ${gtf.name}">' > TMP/cds.hdr


${gtf.name.endsWith(".gz")?"gunzip -c":"cat"} "${gtf}" |\\
	awk -F '\t' '(\$3=="gene") {printf("%s\t%d\t%s\\n",\$1,int(\$4)-1,int(\$5))}' |\\
	sort -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n -T TMP |\\
	bedtools merge   > TMP/gene.bed

echo '##INFO=<ID=in_gene,Number=0,Type=Flag,Description="In gene of ${gtf.name}">' > TMP/gene.hdr


${gtf.name.endsWith(".gz")?"gunzip -c":"cat"} "${gtf}" |\\
	awk -F '\t' '(\$3=="transcript") {printf("%s\t%d\t%s\\n",\$1,int(\$4)-1,int(\$5))}' |\\
	sort -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n -T TMP |\\
	bedtools merge   > TMP/transcript.bed

echo '##INFO=<ID=in_transcript,Number=0,Type=Flag,Description="In transcripts of ${gtf.name}">' > TMP/transcript.hdr

${gtf.name.endsWith(".gz")?"gunzip -c":"cat"} "${gtf}" |\\
	grep -F -w protein_coding |\\
	awk -F '\t' '(\$3=="gene") {printf("%s\t%d\t%s\\n",\$1,int(\$4)-1,int(\$5))}' |\\
	sort -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n -T TMP |\\
	bedtools merge > TMP/protein_coding.bed

echo '##INFO=<ID=in_protein_coding,Number=0,Type=Flag,Description="In protein coding of ${gtf.name}">' > TMP/protein_coding.hdr


bedtools subtract \\
	-a TMP/transcript.bed \\
	-b TMP/exon.bed |\\
	cut -f1,2,3 |\\
	sort -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n -T TMP |\\
	bedtools merge  > TMP/intron.bed

echo '##INFO=<ID=in_intron,Number=0,Type=Flag,Description="In introns of ${gtf.name}">' > TMP/intron.hdr

${gtf.name.endsWith(".gz")?"gunzip -c":"cat"} "${gtf}" |\\
	awk -F '\t' '(\$3=="gene") {S=\$7; B=int(\$4)-1;E=int(\$5);if(S=="+") {E=B-1 ; B = B-${xstream}; if(B<0)B=0; if(E<0) E=0;} else {B=E;E=E+${xstream};} printf("%s\t%d\t%s\\n",\$1,B,E)}' |\\
	sort -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n -T TMP |\\
	bedtools merge > TMP/upstream.bed

echo '##INFO=<ID=in_upstream,Number=0,Type=Flag,Description="Upstream ${xstream}bp of genes in ${gtf.name}">' > TMP/upstream.hdr


${gtf.name.endsWith(".gz")?"gunzip -c":"cat"} "${gtf}" |\\
	awk -F '\t' '(\$3=="gene") {S=\$7; B=int(\$4)-1;E=int(\$5);if(S=="-") {E=B-1 ; B = B-${xstream}; if(B<0)B=0; if(E<0) E=0;} else {B=E;E=E+${xstream};} printf("%s\t%d\t%s\\n",\$1,B,E)}' |\\
	sort -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n -T TMP |\\
	bedtools merge  > TMP/downstream.bed

echo '##INFO=<ID=in_downstream,Number=0,Type=Flag,Description="Downstream ${xstream}bp of genes in ${gtf.name}">' > TMP/downstream.hdr


if ${genes_of_interest?true:false}
then

	cat "${genes_of_interest}" |\\
		grep -vE '^#' |\\
		grep -vE '^\$' |\\
		sort -S ${task.memory.kilo} -t '\t' -k1,1 -T TMP > TMP/jeter.a

	${gtf.name.endsWith(".gz")?"gunzip -c":"cat"} "${gtf}" |\\
		awk -F '\t' '(\$3=="gene")' |\\
		java ${jvm} -jar \${HOME}/jvarkit.jar gtf2bed -c 'gene_name,gtf.strand' |\\
		sort -S ${task.memory.kilo} -t '\t' -k4,4 -T TMP > TMP/jeter.b
	

	join -t '\t' -1 1 -2 4 -o '2.1,2.2,2.3,2.4' TMP/jeter.a TMP/jeter.b | uniq > TMP/roi.bed

	# save roi of gene_name for any future use
	mv TMP/roi.bed ${prefix}.roi.bed

    echo '##INFO=<ID=in_roi,Number=0,Type=Flag,Description="in ROI ${genes_of_interest.name} of ${gtf.name}">' > TMP/roi.hdr

	join -t '\t' -1 1 -2 4 -o '2.1,2.2,2.3,2.5' TMP/jeter.a TMP/jeter.b |\\
		awk -F '\t' '{S=\$4; B=int(\$2)-1;E=int(\$3);if(S=="+") {E=B-1 ; B = B-${xstream}; if(B<0)B=0; if(E<0) E=0;} else {B=E;E=E+${xstream};} printf("%s\t%d\t%s\\n",\$1,B,E)}' |\\
		sort -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n -T TMP |\\
		bedtools merge  > TMP/upstream_roi.bed

	echo '##INFO=<ID=in_roi_upstream,Number=0,Type=Flag,Description="in ROI ${genes_of_interest.name} Upstream ${xstream}bp of genes in ${gtf.name}">' > TMP/roi_upstream.hdr

	join -t '\t' -1 1 -2 4 -o '2.1,2.2,2.3,2.5' TMP/jeter.a TMP/jeter.b |\\
		awk -F '\t' '{S=\$4; B=int(\$2)-1;E=int(\$3);if(S=="-") {E=B-1 ; B = B-${xstream}; if(B<0)B=0; if(E<0) E=0;} else {B=E;E=E+${xstream};} printf("%s\t%d\t%s\\n",\$1,B,E)}' |\\
		sort -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n -T TMP |\\
		bedtools merge  > TMP/downstream_roi.bed

	echo '##INFO=<ID=in_roi_downstream,Number=0,Type=Flag,Description="in ROI ${genes_of_interest.name} Downstream ${xstream}bp of genes in ${gtf.name}">' > TMP/roi_downstream.hdr

fi

bcftools view --threads "${task.cpus}" -O u -o TMP/jeter.bcf "${vcf}"

for F in exon cds gene transcript protein_coding intron upstream downstream roi roi_upstream roi_downstream
do
    if test -f "TMP/\${F}.bed"
    then
		cut -f1,2,3 "TMP/\${F}.bed"  |\\
	        sort -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n -T TMP |\\
			bedtools merge |\\
			sed 's/\$/\t1/' |\\
			bgzip > TMP/jeter.bed.gz
        tabix -f -p bed TMP/jeter.bed.gz

        cat "TMP/\${F}.hdr"  1>&2

        bcftools annotate \\
            ${args1} \\
            --threads "${task.cpus}" \\
            --header-lines "TMP/\${F}.hdr" \\
            --columns "CHROM,POS,END,in_\${F}" \\
            --annotations TMP/jeter.bed.gz \\
            -O u \\
            -o TMP/jeter2.bcf \\
            TMP/jeter.bcf
        
        mv -v TMP/jeter2.bcf TMP/jeter.bcf
    fi
done

bcftools view ${args2}  -O z -o TMP/${prefix}.vcf.gz TMP/jeter.bcf
mv -v TMP/${prefix}.vcf.gz ./

cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
stub:
    def prefix = task.ext.prefix?:"${meta.id}.annot"
"""
touch versions.yml "${prefix}.vcf.gz" "${prefix}.roi.bed"
"""
}