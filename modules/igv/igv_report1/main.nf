/*

Copyright (c) 2025 Pierre Lindenbaum

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
include {verify       } from '../../utils/functions.nf'
include {isBlank      } from '../../utils/functions.nf'

process IGV_REPORT {
tag "${meta.id}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/igv-reports.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta4),path(cytoband)
    tuple val(meta5),path(refgene),path(refgene_tbi)
    tuple val(meta6),path(sample_list)
    tuple val(meta7),path(pedigree)
    tuple val(meta8),path(vcf),path(vcf_idx)
    tuple val(meta ),path(bams),path(bam_idxs)
output:
    tuple val(meta),path("*.html"),emit:html
    tuple val(meta),path("*.index"),emit:index
    path "versions.yml" , emit: versions
script:
    def contig = "${meta.contig?:""}"
    verify(!isBlank(contig),"${task.process} meta.contig is blank")


    def start = "${meta.start?:""}"
    verify(!isBlank(start),"${task.process} meta.start is blank")
    start = (start as int)
    
    def end = "${meta.end?:""}"
    verify(!isBlank(end),"${task.process} meta.end is blank")
    end = (end as int)

    def ref = "${meta.ref?:""}"
    def flanking = (meta.flanking?:(task.ext.flanking?:100)) as int
    verify(flanking>=0 && flanking < 10000,"${task.process} extend out of bound")

    def info_columns= (meta.info?:(task.ext.info?:""))


    def has_vcf = (vcf?true:false)
    if(has_vcf) {
        verify(!isBlank(ref),"${task.process} ref is blank")
        }

    def jvm = task.ext.jvm?:" -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP"
    def title = "${contig}:${start}${start==end?"":"-${end}"}"
    def prefix = task.ext.prefix?:"${meta.id}.igv"
    
"""
hostname 1>&2
mkdir -p  TMP
set -x

# create list of bams
cat << EOF > TMP/jeter.list
${bams.collect{it.name}.join("\n")}
EOF

test -s  TMP/jeter.list

if ${sample_list?true:false}
then
    # keep track of the order of the bams
    samtools samples < TMP/jeter.list |\\
        awk -F '\t' '{printf("%d\t%s\\n",NR,\$0);}' |\\
        sort -T TMP -t '\t' -k2,2 > TMP/jeter.a

    sort -T TMP -t '\t' -k1,1  ${sample_list} > TMP/jeter.b

    # join and restore original order
    join -t '\t' -1 2 -2 1 -o '1.1,1.3' TMP/jeter.a TMP/jeter.b |\\
        sort -t '\t' -k1,1n |\\
        cut -f2 > TMP/jeter.list

    test -s  TMP/jeter.list
fi    

if ${has_vcf}
then

    bcftools view \\
        -i 'REF="${ref}" && POS=${start} && CHROM="${contig}"' \\
        -O z \\
        -o TMP/jeter.vcf.gz \\
        "${vcf}" '${contig}:${start}-${end}'
    
    bcftools index -t TMP/jeter.vcf.gz

    echo "<div id='x1'>" >> TMP/footer.xml
    bcftools view TMP/jeter.vcf.gz |\\
        jvarkit ${jvm} vcf2table \\
            --format html --no-html-header \\
            ${pedigree?"--pedigree ${pedigree}":""} >> TMP/footer.xml
    echo "</div>" >> TMP/footer.xml

fi


# create BED
echo "${contig}\t${(start as int)-1}\t${end}" > TMP/jeter.bed

create_report \\
    TMP/jeter.bed  \\
    ${fasta} \\
    ${cytoband?"--ideogram \"${cytoband}\"":""} \\
    --flanking ${flanking} \\
    ${isBlank(info_columns)?"":"--info-columns \"${info_columns}\""} \\
    --tracks  \\
        ${refgene?"${refgene}":""} \\
        ${has_vcf?"TMP/jeter.vcf.gz":""} \\
        \$(cat TMP/jeter.list) \\
    --output TMP/jeter.html


cat << EOF >  TMP/header.xml
<div>${contig}:${start}${end==start?"":"-"+end}</div>
EOF


echo "<div>" > TMP/footer2.xml
if test -f TMP/footer.xml
then
    cat TMP/footer.xml >> TMP/footer2.xml
fi
mv  TMP/footer2.xml  TMP/footer.xml
cat << EOF >>  TMP/footer.xml
<div>
BAMS: ${bams.collect{"<code>${it.name}</code>"}.join(" , ")}
Generated on \$(date) by \${USER}.
</div>
EOF
echo "</div>" >> TMP/footer.xml

xmllint --noout TMP/footer.xml

cp "${moduleDir}/igv_report.xsl" TMP/jeter.xsl
xmllint --path TMP --xinclude --output TMP/jeter2.xsl  TMP/jeter.xsl
mv  TMP/jeter2.xsl  TMP/jeter.xsl

xsltproc \\
    --xincludestyle \\
    --stringparam title "${contig}:${start}${start==end?"":"-${end}"}" \\
    --html  \\
    TMP/jeter.xsl \\
    TMP/jeter.html > TMP/jeter2.html
mv TMP/jeter2.html TMP/jeter.html


mv -v "TMP/jeter.html" ./${prefix}.html




cat << EOF > ./${prefix}.index
<tr>
    <td>${meta1.ucsc_name?:fasta.name}</th>
    <td><a href="${prefix}.html"><code>${contig}:${start}${end==start?"":"-"+end}</code></a></td>
</tr>
EOF


cat << EOF > versions.yml
    bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
    igvreport: todo
EOF
"""

stub:
    def prefix = task.ext.prefix?:"${meta.id}.igv"
"""
touch versions.yml ${prefix}.html ${prefix}.index
"""
}
