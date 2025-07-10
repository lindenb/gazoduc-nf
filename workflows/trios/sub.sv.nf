
include {TRUVARI                                  } from '../../subworkflows/truvari/main.nf'
include {ANNOTATE_SV                              } from '../../subworkflows/annotation/sv/main.nf'
include {SAMTOOLS_DEPTH_PLOT_COVERAGE             } from  '../../modules/samtools/plotcoverage/main.nf'
include {DELLY                                    } from '../../subworkflows/delly2/main.nf'

workflow WORKFLOW_SV {
take:
    meta
    fasta
    fai
    dict
    bed
    pedigree
    gff3
    gtf
    triosbams_ch
    vcf
main:
    versions = Channel.empty()


    vcf=Channel.empty()

    ch1  = triosbams_ch
        .filter{it[3].equals(fasta[1])}
        .flatMap{[
            [ [id:it[0],status:"case"], it[6], it[7]],
            [ [id:it[1],status:"control"], it[8], it[8]],
            [ [id:it[2],status:"control"], it[9], it[10]]
        ]}
    
    DELLY(meta,fasta,fai,dict,ch1)

    /*
    TRUVARI(
	        meta,
            fasta,
            fai,
            dict,
            vcf
        )
    ANNOTATE_SV(
            meta,
            fasta,
            fai,
            dict,
            gtf,
            TRUVARI.out.vcf
        )

    FILTER_SV(fasta,fai,dict,
        pedigree,
        ANNOTATE_SV.out.vcf
        )

    ch1 = FILTER_SV.out.bed
        .map{it[1]}
        .splitCsv(sep:'\t',header:false)
        .combine(triosbams_ch)
        .filter{it[5].equals(it[6])}
        .map{it.remove(5);return it;}//remove returns deleted item
        .map{[
                [
                id: [it[0], it[1], it[2],it[4],it[3],it[5]].join("_"),
                contig:it[0],
                start:it[1],
                end:it[2],
                svtype:it[3],
                svlen:it[4],
                child:it[5],
                father:it[6],
                mother:it[7],
                title: "SVLEN:"+it[4]+" SVTYPE:"+it[3]+" "+it[5]
                ],
            it[8],//fasta
            it[9],//fai
            it[10],//dict
            [it[11],it[12],it[13],it[14],it[15],it[16]]// bam and bai for child/father/mother
            ]
            }
        .filter{(it[0].svlen as int) > 10000}
        .take(50)
        .view()
    
    SAMTOOLS_DEPTH_PLOT_COVERAGE(gtf,ch1)
    INDEX_FOR_PDFS(
            SAMTOOLS_DEPTH_PLOT_COVERAGE.out.pdf
            .map{it[1].name}
            .collect()
            .map{[[id:"cnv"],it]}
        )
    */
emit:
    versions
}



process FILTER_SV {
tag "${meta.id?:""}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta4),path(pedigree)
    tuple val(meta),path(vcf),path(idx)
output:
    tuple val(meta),path("*.bcf"),path("*.csi"),emit:vcf
    tuple val(meta),path("*.bed"),emit:bed
    path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:vcf.baseName+".filter"
    def minLen = task.ext.min_length?:50
    def maxLen = task.ext.max_length?:100000
"""
set -o pipefail
mkdir -p TMP

# convert pedigree if no 6th column
awk '{S=\$6 ; if(NF==5 || S=="") { if(\$3!="0" && \$4!="0") {S="case";} else {S="control"} }  printf("%s\t%s\t%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,\$4,\$5,S);}' ${pedigree} > TMP/pedigree.tsv

## all other samples are controls
comm -13 \\
	<(cut -f 2 TMP/pedigree.tsv  | sort | uniq) \\
	<(bcftools query -l '${vcf}'| sort | uniq) |\\
	awk '{printf("if(!acceptControl(variant,\\"%s\\")) return false;\\n",\$1);}' >> TMP/custom.m4

awk -F '\t' '(\$6=="control" || \$6=="unaffected") {printf("if(!acceptControl(variant,\\"%s\\")) return false;\\n",\$2);}'  TMP/pedigree.tsv >> TMP/custom.m4

awk -F '\t' '(\$6=="case" || \$6=="affected") {printf("if(acceptTrio(variant,\\"%s\\",\\"%s\\",\\"%s\\")) children.add(\\"%s\\");\\n",\$2,\$3,\$4,\$2);}' TMP/pedigree.tsv  >> TMP/custom.m4


m4 -P -D__MIN_LEN__=${minLen} -D__MAX_LEN__=${maxLen} -I TMP < "${moduleDir}/select.sv.m4"  > TMP/jeter.code

touch TMP/jeter.bed

bcftools view ${vcf} |\\
    jvarkit -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP vcffilterjdk \\
        --body  -f  TMP/jeter.code |\\
    bcftools view --write-index -O b -o TMP/jeter.bcf


LC_ALL=C sort -t '\t' -T TMP -k1,1 -k2,2n TMP/jeter.bed | uniq > TMP/jeter2.bed


mv TMP/jeter.bcf ${prefix}.bcf
mv TMP/jeter.bcf.csi ${prefix}.bcf.csi
mv TMP/jeter2.bed ${prefix}.bed

cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
}

process INDEX_FOR_PDFS {
tag "${meta.id?:""}"
label "process_single"
afterScript "rm -rf TMP"
input:
    tuple val(meta),val(names)
output:
    tuple val(meta),path("index.html"),emit:html
script:
"""
cat << EOF > index.html
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
<meta http-equiv="author" content="Pierre Lindenbaum Phd ">
<title>CNV De NOVO</title>
<script>
var files=[${names.collect{"\"${it}\""}.join(",")}];
var page =0;

function goTo(dx) {
    if(files.length==0) return;
    page = page+dx;
    if(page<0) page = files.length-1;
    if(page>=files.length) page=0;
    document.getElementById("id1").src = files[page];
    document.getElementById("h1").textContent = files[page]+ " ("+(page+1)+"/"+files.length+")";
    }

 

window.addEventListener('load', (event) => {
  var frame = document.getElementById("id1");
  frame.style.height=(frame.contentWindow.document.body.scrollHeight+20)+'px';
  goTo(0);
});

</script>
</head>
<body>
<div>
    <button onclick="goTo(-1);">PREV</button>
    <span id="h1"></span>
    <button onclick="goTo(1);">NEXT</button>
</div>

<iframe id="id1" style="height:800px;width:100%;" src="blank_">
</iframe>

<div>
    <button onclick="goTo(-1);">PREV</button>
    <span>navigation</span>
    <button onclick="goTo(1);">NEXT</button>
</div>


</body>
</html>
EOF
"""
}