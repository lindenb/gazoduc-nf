
/*

THIS FILE WAS GENERATED DO NOT EDIT !!!

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

include { flatMapByIndex                      } from '../../../modules/utils/functions.nf'
include { makeKey                             } from '../../../modules/utils/functions.nf'
include { ALPHAGENOME_BATCH_VARIANT_SCORING   } from '../../../modules/alphagenome/batch_variant_scoring'


workflow APHAGENOME_VARIANT_SCORING {
take:  
    metadata
    fasta
    fai
    dict
    vcfs
main:
    versions =  Channel.empty()
    SPLIT_VCF(vcfs)
    versions = versions.mix(SPLIT_VCF.out.versions)
    
    ALPHAGENOME_BATCH_VARIANT_SCORING(
        SPLIT_VCF.out.tsv
            .flatMap{row->flatMapByIndex(row,1)}
            .map{meta,tsv->[meta.plus(id:tsv.toRealPath().toString().md5(),old_id:meta.id),tsv]}
        )
     versions = versions.mix(ALPHAGENOME_BATCH_VARIANT_SCORING.out.versions)
    /*
    track_name_ch = ALPHAGENOME_BATCH_VARIANT_SCORING.out.tsv
        .map{_meta,tsv->tsv}
        .splitCsv(header:true,sep:'\t')
        .map{row->row.track_name}
        .filter{track_name->track_name=="donor" || track_name=="acceptor"}
        .unique()
        .view()
    */
    track_name_ch = Channel.of("donor","acceptor")

    score_ch = Channel.of(
        ["raw_score",24],
        ["quantile_score",25]
        ).map{k,v->[score_name:k,score_column:v]}

    TOXML(
        dict,
        ALPHAGENOME_BATCH_VARIANT_SCORING.out.tsv
            .combine(score_ch)
            .map{meta1,tsv,meta2->[meta1.plus(meta2),tsv]}
            .combine(track_name_ch)
            .map{meta1,tsv,track_name_s->[meta1.plus(track_name:track_name_s),tsv]}
            .map{meta,tsv->[[track_name:meta.track_name, score_name: meta.score_name, score_column: meta.score_column],tsv]}
            .map{meta,tsv->[meta.plus(id:makeKey(meta)),tsv]}
            .groupTuple()
            .map{meta,files->[meta,files.sort()]}
        )
    versions = versions.mix(TOXML.out.versions)

emit:
    versions
}

process TOXML {
label "process_single"
tag "${meta.id} ${meta.track_name} ${meta.score_name}"
label "process_single"
afterScript "rm -rf TMP"
conda "../../gazoduc-nf/conda/bioinfo.01.yml"
input:
    tuple val(meta2),path(dict)
    tuple val(meta ),path(tsvs)
output:
    tuple val(meta ),path("*.xml.gz"),emit:xml
    path("versions.yml"),emit:versions
script:
    def score_column = meta.score_column
    def track_name = meta.track_name
    def jvm =  task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData"
    def jvarkit = task.ext.jvarkit?:"java ${jvm} -jar \${HOME}/jvarkit.jar"// TODO update when conda released
    def prefix = task.ext.prefix?:"${meta.id}.bed2xml"
"""
mkdir -p TMP
gunzip -c ${tsvs} |\\
    awk -F '\t' '(\$11=="${track_name}" && \$1!="variant_id") {split(\$1,a,/[:]/); printf("%s\t%d\t%s\t%s\\n",a[1],int(a[2])-1,a[2],\$${score_column});}' |\\
     ${jvarkit} bed2xml  \\
        --regex "(chr)?[0-9XY]+" \\
        --type bedGraph \\
        --columns score \\
        -R ${dict} |gzip --best > TMP/${prefix}.xml.gz

mv TMP/${prefix}.xml.gz ./

cat << EOF > versions.yml
${task.process}:
	jvarkit: "\$(jvarkit --version)"
    xmllint: \$(xmllint --version 2>&1 |awk '(NR==1) {print \$NF;}')
EOF
"""

stub:
    def prefix = task.ext.prefix?:"${meta.id}.bed2xml"
"""
touch versions.yml ${prefix}.xml.gz
"""
}

process SPLIT_VCF {
tag "${meta.id}"
label "process_single"
afterScript "rm -rf TMP"
conda "../../gazoduc-nf/conda/bioinfo.01.yml"
input:
        tuple val(meta),path(vcf)
output:
        tuple val(meta),path("*.tsv", arity: '1..*'),emit:tsv
        path("versions.yml"),emit:versions
script:
        def nlines = task.ext.nlines?:"100"
        def prefix = task.ext.prefix?:"${meta.id}.split."

//      java -jar \${HOME}/jvarkit.jar vcfgnomad  --gnomad "${params.gnomad}" --fields "AF_nfe"  --max-af '0.01' |\\
//      bcftools view --apply-filters '.,PASS' -O u |\\

"""
mkdir -p TMP
bcftools norm  -m -any -O u  "${vcf}" |\\
        bcftools query -f '%VKX\t%CHROM\t%POS\t%REF\t%ALT\\n' |\\
        awk -F '\t' '((\$4 ~ /^[ACGTN]+\$/ && \$5 ~ /^[ACGTN]+\$/) && (\$2 ~ /^(chr)?[0-9XY]*\$/))' |\\
        sort -T TMP | uniq |\
        split --additional-suffix=.tsv --lines "${nlines}" - TMP/${prefix}

find TMP/ -type f -name "*.tsv" | while read F
do
        echo -e 'variant_id\tCHROM\tPOS\tREF\tALT' > TMP.txt
        cat "\$F" >> TMP.txt
        mv TMP.txt "\$F"
done

mv TMP/*.tsv ./ || true

touch versions.yml
"""

stub:
"""
touch versions.yml
"""
}
