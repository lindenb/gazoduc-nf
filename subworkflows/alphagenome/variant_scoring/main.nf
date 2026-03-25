
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
include { ZIP                                 } from '../../../modules/utils/zip'
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
    SPLIT_VCF(
        [[id:"blacklist"],file("${moduleDir}/blacklisted.tsv")],
        vcfs
        )
    versions = versions.mix(SPLIT_VCF.out.versions)
    
    ALPHAGENOME_BATCH_VARIANT_SCORING(
        SPLIT_VCF.out.tsv
            .flatMap{row->flatMapByIndex(row,1)}
            .map{meta,tsv->[meta.plus(id:makeKey(tsv),old_id:meta.id),tsv]}
        )
     versions = versions.mix(ALPHAGENOME_BATCH_VARIANT_SCORING.out.versions)

    output_type_ch = Channel.of("ATAC" ,"CAGE" ,"CHIP_HISTONE" ,"CHIP_TF" ,
                "CONTACT_MAPS" ,"DNASE" ,"PROCAP" ,"RNA_SEQ" ,"SPLICE_JUNCTIONS" ,
                "SPLICE_SITES" ,"SPLICE_SITE_USAGE"
                );


    score_ch = Channel.of(
        ["raw_score",24],
        ["quantile_score",25]
        ).map{k,v->[score_name:k,score_column:v]}
   

    by_output_type_ch = ALPHAGENOME_BATCH_VARIANT_SCORING.out.tsv
            .combine(score_ch)
            .map{meta1,tsv,meta2->[meta1.plus(meta2),tsv]}
            .combine(output_type_ch)
            .map{meta1,tsv,output_type_s->[meta1.plus(output_type:output_type_s),tsv]}
            .map{meta,tsv->[[output_type:meta.output_type, score_name: meta.score_name, score_column: meta.score_column],tsv]}
            .map{meta,tsv->[meta.plus(id:makeKey(meta)),tsv]}
            .groupTuple()
            .map{meta,files->[meta,files.sort()]}
           

    PLOT(
        dict,
        by_output_type_ch
        )

    MAKE_ANNOT_FILES(
        by_output_type_ch.filter{meta,_files->meta.score_name=="raw_score"}
        )

   

    ANNOT_VCF(
        vcfs,
         MAKE_ANNOT_FILES.out.tabix
        .map{meta,tsv,idx,hdr,cols->[[id:"annot"],meta,tsv,idx,hdr,cols]}
        .groupTuple()
        .view()
        )

    ZIP(
       PLOT.out.png.map{meta,png->png}.collect().map{f->[[id:"alphagenome"],f.sort()]} 
    )    
    // versions = versions.mix(TOXML.out.versions)

emit:
    versions
}

process PLOT {
label "process_single"
tag "${meta.id} ${meta.output_type} ${meta.score_name}"
label "process_single"
afterScript "rm -rf TMP"
conda "../../gazoduc-nf/conda/bioinfo.02.yml"
input:
    tuple val(meta2),path(dict)
    tuple val(meta ),path(tsvs)
output:
    tuple val(meta ),path("*.png"),emit:png
    path("versions.yml"),emit:versions
script:
    def score_column = meta.score_column
    def output_type = meta.output_type
    def prefix = task.ext.prefix?:"${meta.id}.${meta.score_name}.${output_type}"
    def jvm =  task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData"
    def jvarkit = task.ext.jvarkit?:"java ${jvm} -jar \${HOME}/jvarkit.jar"// TODO update when conda released
    def plot_size = 400
"""
mkdir -p TMP
set -x
gunzip -c ${tsvs} |\\
    awk -F '\t' '(\$9=="${output_type}" && \$1!="variant_id") {split(\$1,a,/[:]/); printf("%s\t%d\t%s\t%s\\n",a[1],int(a[2])-1,a[2],\$${score_column});}' |\\
    ${jvarkit} addlinearindextobed -R ${dict} --ignore --regex '(chr)?[0-9XY]+' |\\
    cut -f1,5 |\\
    sort -S '${task.memory.kilo}' -t '\t' -T TMP -k1,1n -k2,2g  > TMP/jeter.tsv

 ${jvarkit} dict2bed --no-header ${dict}  |\\
    ${jvarkit} addlinearindextobed -R ${dict} --ignore --regex '(chr)?[0-9XY]+' |\\
    awk '{printf("abline(v = %s.0)\\n",\$1);}' > TMP/contigs.R



cat << '__EOF__' | R --no-save
d <- read.table("TMP/jeter.tsv", header=FALSE, stringsAsFactors=FALSE,col.names = c("POS", "SCORE"), colClasses = c("numeric","numeric"))

head(d)

png("TMP/jeter.png", width = ${(plot_size as int)*4}, height = ${plot_size}, unit = "px")
plot(
    x = d\$POS,
    y = d\$SCORE,
    main="${output_type} ${meta.score_name}",
    xlab="genomic position",
    ylab="${meta.score_name}"
    )
source("TMP/contigs.R",local=TRUE)
dev.off()
__EOF__

mv TMP/jeter.png "${prefix}.png"

cat << EOF > versions.yml
${task.process}:
EOF
"""

stub:
    def prefix = task.ext.prefix?:"${meta.id}"
"""
touch versions.yml ${prefix}.png
"""
}

process MAKE_ANNOT_FILES {
label "process_single"
tag "${meta.id} ${meta.output_type} ${meta.score_name}"
label "process_single"
afterScript "rm -rf TMP"
conda "../../gazoduc-nf/conda/bioinfo.01.yml"
input:
    tuple val(meta ),path(tsvs)
output:
    tuple val(meta ),path("*.tsv.gz"),path("*.tsv.gz.tbi"),path("*.hdr"),path("*.columns"),emit:tabix
    path("versions.yml"),emit:versions
script:
    def score_column = meta.score_column
    def prefix = task.ext.prefix?:"${meta.id}.${meta.score_name}.${meta.output_type}"
"""
mkdir -p TMP
set -x
gunzip -c ${tsvs} |\\
    awk -F '\t' '(\$9=="${meta.output_type}" && \$1!="variant_id") {printf("%s\t%s\\n",\$1,\$${score_column});}' |\\
    LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2gr |\\
    awk -F '\t' 'BEGIN{PREV=""} {if(\$1==PREV) next; print;PREV=\$1;}' |\\
    awk -F '\t' '{split(\$1,a,/[:>]/); printf("%s\t%s\t%s\t%s\t%s\\n",a[1],int(a[2]),a[3],a[4],\$2);}' |\\
    sort -T TMP -t '\t' -k1,1 -k2,2n -k3,3 -k4,4 |\\
    bgzip > TMP/jeter.tsv.gz

tabix -f -s 1 -b 2 -e 2 TMP/jeter.tsv.gz


echo '##INFO=<ID=AG_${meta.output_type},Number=1,Type=Float,Description="ALPHA GENOME max ${meta.score_name}/${meta.output_type}">' > ${prefix}.hdr
echo 'CHROM,POS,REF,ALT,INFO/AG_${meta.output_type}' > ${prefix}.columns

mv TMP/jeter.tsv.gz ${prefix}.tsv.gz
mv TMP/jeter.tsv.gz.tbi ${prefix}.tsv.gz.tbi
touch versions.yml
"""
stub:
def prefix = task.ext.prefix?:"${meta.id}.${meta.score_name}.${meta.output_type}"
"""
touch versions.yml  ${prefix}.tsv.gz 
"""
}


process TOXML {
label "process_single"
tag "${meta.id} ${meta.output_type} ${meta.score_name}"
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
    def output_type = meta.output_type
    def jvm =  task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData"
    def jvarkit = task.ext.jvarkit?:"java ${jvm} -jar \${HOME}/jvarkit.jar"// TODO update when conda released
    def prefix = task.ext.prefix?:"${meta.id}.bed2xml"
"""
mkdir -p TMP
gunzip -c ${tsvs} |\\
    awk -F '\t' '(\$9=="${output_type}" && \$1!="variant_id") {split(\$1,a,/[:]/); printf("%s\t%d\t%s\t%s\\n",a[1],int(a[2])-1,a[2],\$${score_column});}' |\\
     ${jvarkit} bed2xml  \\
        --regex "(chr)?[0-9XY]+" \\
        --type bedGraph \\
        --columns score \\
        -R ${dict} |gzip --best > TMP/${prefix}.xml.gz

mv TMP/${prefix}.xml.gz ./

cat << EOF > versions.yml
${task.process}:
	jvarkit: "\$(${jvarkit} --version)"
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
        tuple val(meta1 ),path(blacklisted)
        tuple val(meta ) ,path(vcf)
output:
        tuple val(meta),path("*.tsv", arity: '1..*'),emit:tsv
        path("versions.yml"),emit:versions
script:
        def args1 = task.ext.args1?:""
        def args2 = task.ext.args2?:""
        def nlines = task.ext.nlines?:"100"
        def prefix = task.ext.prefix?:"${meta.id}.split."

//      java -jar \${HOME}/jvarkit.jar vcfgnomad  --gnomad "${params.gnomad}" --fields "AF_nfe"  --max-af '0.01' |\\
//      bcftools view --apply-filters '.,PASS' -O u |\\

"""
mkdir -p TMP

sort -T TMP -S ${task.memory.kilo}  ${blacklisted} > TMP/exclude.txt

bcftools norm  ${args1} -m -any -O u  "${vcf}" |\\
        bcftools query ${args2} -f '%VKX\t%CHROM\t%POS\t%REF\t%ALT\\n' |\\
        awk -F '\t' '((\$4 ~ /^[ACGTN]+\$/ && \$5 ~ /^[ACGTN]+\$/) && (\$2 ~ /^(chr)?[0-9XY]*\$/))' |\\
        sort -T TMP -S ${task.memory.kilo} | uniq |\
        comm  -3 -2 -  TMP/exclude.txt |\
        split --additional-suffix=.tsv --lines "${nlines}" - TMP/${prefix}

find TMP/ -type f -name "*.tsv" | while read F
do
        echo -e 'variant_id\tCHROM\tPOS\tREF\tALT' > TMP/_jeter.txt
        cat "\$F" >> TMP/_jeter.txt
        mv TMP/_jeter.txt "\$F"
done

mv TMP/*.tsv ./ || true

touch versions.yml
"""

stub:
"""
touch versions.yml
"""
}

process ANNOT_VCF {
tag "${meta.id}"
label "process_single"
afterScript "rm -rf TMP"
conda "../../gazoduc-nf/conda/bioinfo.01.yml"
input:
        tuple val(meta ) ,path(vcf)
        tuple val(meta1 ),val(metas),path(tsvs),path(tibis),path(hdrs),path(columns)

output:
        tuple val(meta),path("*.vcf.gz"),emit:vcf
script:
    def prefix="${meta.id}.alphagenome"
"""
mkdir -p TMP
set -x
bcftools view --write-index -O z -o TMP/jeter1.vcf.gz ${vcf}

find -name "*.tsv.gz" | while read F
do

    bcftools annotate \\
        -a \${F} \\
        --columns-file `basename \${F} .tsv.gz`.columns \
        --header-lines `basename \${F} .tsv.gz`.hdr \\
         -O z -o TMP/jeter2.vcf.gz \\
         TMP/jeter1.vcf.gz

    mv TMP/jeter2.vcf.gz TMP/jeter1.vcf.gz 

    bcftools index  -f -t TMP/jeter1.vcf.gz

done

mv TMP/jeter1.vcf.gz ${prefix}.vcf.gz
"""
}