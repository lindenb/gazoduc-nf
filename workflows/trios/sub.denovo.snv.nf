include {BCFTOOLS_CONCAT                           } from '../../modules/bcftools/concat/main.nf'
include {DOWNLOAD_CYTOBAND                        } from '../../modules/ucsc/download.cytobands/main.nf'
include {DOWNLOAD_REFGENE                         } from '../../modules/ucsc/download.refgene/main.nf'


workflow WORKFLOW_DENOVO_SNV {
take:
    meta
    fasta
    fai
    dict
    gff3
    gtf
    triosbams_ch
    pedigree
    vcf
main:
    versions = Channel.empty()
    
    KEEP_DE_NOVO(pedigree,vcf)
    versions =  versions.mix(KEEP_DE_NOVO.out.versions)

    BCFTOOLS_CONCAT(
        KEEP_DE_NOVO.out.vcf
            .map{[it[1],it[2]]}
            .collect()
            .map{[meta,it.flatten()]},
        [[id:"nobed"],[]
        ])
    versions =  versions.mix(KEEP_DE_NOVO.out.versions)
    vcf = BCFTOOLS_CONCAT.out.vcf



    REPORT(pedigree,vcf)
    versions =  versions.mix(REPORT.out.versions)

    ch1 = REPORT.out.bed
        .map{it[1]}
        .splitCsv(sep:'\t',header:false)
        .combine(triosbams_ch)
        .filter{it[3].equals(it[5])}
        .map{it.remove(5);return it;}//remove returns deleted item
        .map{it.add(0,[id:it[0]+"_"+it[1]+"_"+it[2]+"_"+it[3]]);return it;}
        //.view()
    

    DOWNLOAD_REFGENE(fasta,fai,dict)
    versions =  versions.mix(DOWNLOAD_REFGENE.out.versions)
    DOWNLOAD_CYTOBAND(fasta,fai,dict)
    versions =  versions.mix(DOWNLOAD_CYTOBAND.out.versions)

    IGV_REPORTS(
        DOWNLOAD_CYTOBAND.out.output,
        DOWNLOAD_REFGENE.out.output,
        REPORT.out.vcf,
        ch1
    )
    versions =  versions.mix(IGV_REPORTS.out.versions)

    GATHER_IGV_REPORTS(
        IGV_REPORTS.out.html.map{it[1]}.collect().map{[[id:"igvregport"],it]},
        IGV_REPORTS.out.index.map{it[1]}.collect().map{[[id:"igvregport"],it]}
        )

emit:
    vcf = REPORT.out.vcf
    versions
}



process KEEP_DE_NOVO {
    tag "${meta.id}"
    label "process_single"
    conda "${moduleDir}/../../conda/bioinfo.01.yml"
    afterScript "rm -rf TMP"
    input:
        tuple val(met1),path(pedigree)
        tuple val(meta),path(vcf),path(vcfidx)
    output:
        tuple val(meta),path("*.bcf"),path("*.csi"),optional:true,emit:vcf
        path("versions.yml"),emit:versions
    script:
        def prefix = task.ext.prefix?:vcf.baseName.md5().substring(0,7) +".denovoonly"
    """
    mkdir -p TMP
    
    bcftools view ${vcf} |\\
        jvarkit -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP vcffilterjdk \\
            -e 'return variant.hasAttribute(\\"loConfDeNovo\\")|| variant.hasAttribute(\\"hiConfDeNovo\\") || variant.getAttributeAsInt(\\"MERR\\",0)>0;' |\\
            awk '/#CHROM/ {printf("##INFO=<ID=CONTROLS_HAVING_ALT,Number=1,Type=Integer,Description=\\"In the pedigree, number of control or parent that have ALT allele.\\">\\n");} {print}' |\\
            bcftools view  -O z -o TMP/jeter.jeter.vcf.gz
    
    echo "Genotype g=null; int count=0;" > TMP/jeter.code
    awk '(\$6=="control" || (\$3=="0" && \$4=="0")) {printf("%s\\n",\$2);} ' '${pedigree}'  |\\
        sort | uniq |\\
        awk '{printf("g=variant.getGenotype(\\"%s\\"); if(g!=null && g.hasAltAllele()) count++;\\n",\$1);}' >> TMP/jeter.code
        echo 'return new VariantContextBuilder(variant).attribute("CONTROLS_HAVING_ALT",count).make();' >>  TMP/jeter.code
    
    bcftools view TMP/jeter.jeter.vcf.gz |\\
        jvarkit -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP vcffilterjdk \\
            -f  TMP/jeter.code |\\
            bcftools view  -O z -o TMP/jeter.jeter2.vcf.gz
    
    bcftools view --write-index -O b -o TMP/jeter.bcf TMP/jeter.jeter2.vcf.gz

    rm TMP/jeter.jeter.vcf.gz
    rm TMP/jeter.jeter2.vcf.gz

    if test \$(bcftools index -s TMP/jeter.bcf |wc -l) -gt 0
    then
        mv TMP/jeter.bcf ${prefix}.bcf
        mv TMP/jeter.bcf.csi ${prefix}.bcf.csi
    fi

cat << END_VERSIONS > versions.yml
"${task.process}":
    bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
    """
}

process REPORT {
    tag "${meta.id}"
    label "process_single"
    conda "${moduleDir}/../../conda/bioinfo.01.yml"
    afterScript "rm -rf TMP"
    input:
        tuple val(met1),path(pedigree)
        tuple val(meta),path(vcf),path(vcfidx)
    output:
        tuple val(meta),path("*.vcf.gz"),path("*.tbi"),emit:vcf
        tuple val(meta),path("*.table.txt")
        tuple val(meta),path("*.genes.tsv")
        tuple val(meta),path("*.bed"),optional:true,emit:bed //used for igv reports
        path("versions.yml"),emit:versions
    script:
        def prefix = task.ext.prefix?:"snv.denovo"
        def args1 = task.ext.args1?:""
	// do not create IGV report for more than XX variants
	def max_igv_reports= task.ext.max_igv_reports?:1000
    """
    mkdir -p TMP
    set -o pipefail


    bcftools view  ${args1} -Oz -o TMP/jeter.vcf.gz "${vcf}"
    bcftools index -f -t TMP/jeter.vcf.gz 
   
cat << '_EOF_' > TMP/jeter.awk
    {
    for(x=4;x<=5;++x) {
        N=split(\$x,a,/[,]/);
        for(i=1;i<=N;i++) {
            if(a[i]=="." || a[i]=="") continue;
            printf("%s\t%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,a[i],(x==4?"loConfDeNovo":"hiConfDeNovo"));
            }
        }
    }
_EOF_


    bcftools query -f '%CHROM\t%POS0\t%END\t%INFO/loConfDeNovo\t%INFO/hiConfDeNovo\\n' TMP/jeter.vcf.gz |\\
        awk -F '\t' -f TMP/jeter.awk |\\
        LC_ALL=C sort -T TMP -k1,1 -k2,2n |\\
        uniq > TMP/jeter.bed
    
    bcftools view TMP/jeter.vcf.gz |\\
        jvarkit -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP vcf2table \\
            --hide 'HOM_REF,NO_CALL' --pedigree ${pedigree} > TMP/jeter.table.txt 

    bcftools view TMP/jeter.vcf.gz |\\
        jvarkit -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP groupbygene > TMP/jeter.genes.tsv


    mv TMP/jeter.vcf.gz ./${prefix}.vcf.gz
    mv TMP/jeter.vcf.gz.tbi ./${prefix}.vcf.gz.tbi
    mv TMP/jeter.table.txt ./${prefix}.table.txt
    mv TMP/jeter.genes.tsv ./${prefix}.genes.tsv


    if test \$(cat TMP/jeter.bed |wc -l) -lt ${max_igv_reports}
    then
        grep -vFw loConfDeNovo  TMP/jeter.bed  > TMP/jeter2.bed
	mv TMP/jeter2.bed TMP/jeter.bed
    fi


    if test \$(cat TMP/jeter.bed |wc -l) -lt ${max_igv_reports}
    then
        mv TMP/jeter.bed  ./${prefix}.bed
    fi


cat << END_VERSIONS > versions.yml
"${task.process}":
    bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
    """
}



process IGV_REPORTS {
tag "${meta.id?:""}"
label "process_single"
conda "${moduleDir}/../../conda/igv-reports.yml"
afterScript "rm -rf TMP"
input:
        tuple val(meta4),path(cytoband)
        tuple val(meta5),path(refgene),path(refgene_tbi)
        tuple val(meta6),path(vcf),path(vcfidx)
        tuple val(meta),val(contig),val(start0),val(end0),
                val(child),val(quality),val(father),val(mother),
                path(fasta),path(fai),path(dict),
                path(bamC),path(baiC),
                path(bamP),path(baiP),
                path(bamM),path(baiM)
output:
        tuple val(meta),path("*.html"),emit:html
        tuple val(meta),path("*.index"),emit:index
        path("versions.yml"),emit:versions
script:
    def flanking = task.ext.flanking?:100
    def start = ((start0 as int)+1)
    def end = (end0 as int)
    def info_columns= task.ext.info?:"CARDIOPANEL CLINVAR_CLNSIG REVEL BCSQ ANN MERR hiConfDeNovo loConfDeNovo,ALPHAMISSENSE_PATHOGENOCITY ALPHAMISSENSE_CLASS"
    def prefix = contig+  "_" + String.format("%09d",start) + (start==end?"":"_"+String.format("%09d",end)) +"_" + child
"""
hostname 1>&2
mkdir -p TMP

bcftools view \\
    --samples "${child},${father},${mother}" \\
    --trim-unseen-allele \\
    --trim-alt-alleles \\
    --force-samples \\
    --regions "${contig}:${start}-${end}" \\
    -O u \\
    "${vcf}" |\\
    bcftools view -i 'AC[*]>0' -O z -o TMP/jeter.vcf.gz

bcftools index -f -t TMP/jeter.vcf.gz

echo -e "${contig}\t${start0}\t${end0}" > TMP/jeter.bed

##gunzip -c "${refgene}" > TMP/${refgene.baseName}

cat << EOF >  TMP/header.html
<div>
De novo variant at ${contig}:${start}${end==start?"":"-"+end} for <b>child</b>:<i>${child}</i> , <b>father</b>:<i>${father}</i> , <b>mother</b>:<i>${mother}</i>
</div>
EOF


cat << EOF >  TMP/footer.html
<div>
Bam files: <code>${bamC}</code>, <code>${bamP}</code> , <code>${bamM}</code>.<br/>
Generated on \$(date) by \${USER}.
</div>
EOF

create_report TMP/jeter.bed  \\
    ${fasta} \\
    --header TMP/header.html \\
    --header TMP/footer.html \\
	--ideogram "${cytoband}" \\
	--flanking ${flanking} \\
	${info_columns.isEmpty()?"":"--info-columns \"${info_columns}\""} \\
	--tracks  TMP/jeter.vcf.gz ${bamC} ${bamP} ${bamM} ${refgene} \\
	--output TMP/jeter.html

mv -v "TMP/jeter.html" ./${prefix}.html

cat << EOF > TMP/jeter.html
<tr class="${quality}">
    <td><code>${contig}:${start}${end==start?"":"-"+end}</code></td>
    <td>(${fasta.baseName})</th>
    <td><a href="${prefix}.html">${child}</a></td>
    <td>${father}</td>
    <td>${mother}</td>
    <td><span class="${quality}">${quality}</span></td>
</tr>
EOF

mv -v "TMP/jeter.html" ./${prefix}.index

cat << END_VERSIONS > versions.yml
"${task.process}":
    bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
    igvreport: todo
END_VERSIONS
"""
}

process GATHER_IGV_REPORTS {
    tag "${meta.id?:""}"
label "process_single"
conda "${moduleDir}/../../conda/igv-reports.yml"
afterScript "rm -rf TMP"
input:
    tuple val(meta),path("PAGES/*")
    tuple val(meta2),path("INDEX/*")
output:
    tuple val(meta),path("index.html"),emit:html
    tuple val(meta),path("*.zip"),emit:zip
    path("versions.yml"),emit:versions
script:
    def title = task.ext.title?:"De Novo"
    def prefix = task.ext.prefix?:"archive"
"""
cat << EOF > index.html
<html>
<head>
    <title>${title}</title>
<style>
table, th, td {
  border: 1px solid black;
  border-collapse: collapse;
}

tr.hiConfDeNovo {
  background-color:#C7DBF0;
}

tr.loConfDeNovo {
  background-color:yellow;
}

</style>
</head>
<body>
<table>
<thead>
    <tr>
        <th>Position</th>
        <th>build</th>
        <th>Child</th>
        <th>Father</th>
        <th>Mother</th>
        <th>Quality</th>
    </tr>
</thead>
<tbody>
EOF

find INDEX/ -name "*.index" |\
    sort -T . -V |\
    xargs -L10 cat >> index.html

cat << EOF >> index.html
</tbody>
</table>
</body>
</html>
EOF

mkdir -p "${prefix}"
cp index.html "${prefix}/"
cp PAGES/*.html "${prefix}/"

zip -r9 ${prefix}.zip  "${prefix}/"

touch versions.yml
"""
}
