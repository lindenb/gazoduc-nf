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
include {MULTIQC                     } from '../../modules/multiqc'
include {COMPILE_VERSIONS            } from '../../modules/versions/main.nf'
include {runOnComplete               } from '../../modules/utils/functions.nf'
include {isBlank                     } from '../../modules/utils/functions.nf'
include {DRAGEN_SCAN_DIRECTORY       } from '../../modules/dragen/scan.directory'
include {PREPARE_ONE_REFERENCE       } from '../../subworkflows/samtools/prepare.one.ref'
include {DOWNLOAD_GNOMAD_SV          } from '../../modules/gnomad_sv/download.vcf/main.nf'

workflow {
    versions = Channel.empty()
    multiqc = Channel.empty()
    metadata = [id: "dragen"]

    if(isBlank(params.base)) {
        log.warn("undefined --base");
        exit -1
        }
    if(isBlank(params.fasta)) {
        log.warn("undefined --fasta");
        exit -1
        }
    DRAGEN_SCAN_DIRECTORY([metadata,file(params.base)])
    versions = versions.mix(DRAGEN_SCAN_DIRECTORY.out.versions)

    dispatch = DRAGEN_SCAN_DIRECTORY.out.json
        .map{_meta,json->json}.splitJson().branch{v->
        ploidy: v.type=="ploidy"
        snv: v.type=="vcf"
        sv : v.type=="sv"
        other:true
        }

    PLOIDY(dispatch.ploidy.map{[[id:it.id],file(it.vcf)]}) 
    versions = versions.mix(PLOIDY.out.versions)



  /***************************************************
   *
   *  PREPARE FASTA REFERENCE
   *
   */
	PREPARE_ONE_REFERENCE(
		metadata,
		Channel.fromPath(params.fasta).map{[[id:it.baseName],it]}
		)
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)


    PREPARE_BED(
        PREPARE_ONE_REFERENCE.out.fasta,
        PREPARE_ONE_REFERENCE.out.fai,
        PREPARE_ONE_REFERENCE.out.dict,
        [[id:"gtf"],file(params.gtf)],
        [[id:"genes"],file(params.genes)]
        )
    versions = versions.mix(PREPARE_BED.out.versions)


    ANNOT_SNV(
        PREPARE_ONE_REFERENCE.out.fasta,
        PREPARE_ONE_REFERENCE.out.fai,
        PREPARE_ONE_REFERENCE.out.dict,
        PREPARE_BED.out.bed,
        dispatch.snv
		.map{[[id:it.id],file(it.vcf),file(it.tbi)]}
		.filter{_meta,f,_t->f.parent!=null  && f.parent.name!="VCF"}
        )
     versions = versions.mix(ANNOT_SNV.out.versions)

    

    DOWNLOAD_GNOMAD_SV(
        PREPARE_ONE_REFERENCE.out.dict
        )
     versions = versions.mix(DOWNLOAD_GNOMAD_SV.out.versions)

    ANNOT_SV(
        PREPARE_ONE_REFERENCE.out.fasta,
        PREPARE_ONE_REFERENCE.out.fai,
        PREPARE_ONE_REFERENCE.out.dict,
        PREPARE_BED.out.tabix,
        DOWNLOAD_GNOMAD_SV.out.vcf,
        dispatch.sv.map{[[id:it.id],file(it.vcf),file(it.tbi)]}
        )
    versions = versions.mix(ANNOT_SV.out.versions)

    MERGE_VARIANTS(
        ANNOT_SNV.out.vcf.map{meta,vcf,tbi->vcf}.collect().map{[[id:"snv"],it.sort()]}.mix(
                ANNOT_SV.out.vcf.map{meta,vcf,tbi->vcf}.collect().map{[[id:"sv"],it.sort()]}
                )
    )

}
process PLOIDY {
label "process_single"
label "array100"
tag "${meta.id}"
afterScript "rm -rf TMP"
conda   "${moduleDir}/../../conda/bioinfo.01.yml"
input:
    tuple val(meta),path(vcf)
output:
    tuple val(meta),path("*.tsv"),optional:true,emit:tsv
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:vcf.name
"""
mkdir -p TMP

bcftools view --header-only "${vcf}" |\\
    awk -F '=' '
        BEGIN {P="";} 
        /^##estimatedSexKaryotype=/ {P=\$2;}
        END {if(P!="XX" && P="XY") {printf("${meta.id}\t%s\\n",P);}}
        ' > TMP/jeter.txt

if test -s TMP/jeter.txt
then
    mv TMP/jeter.txt "${prefix}.ploidy.tsv"
fi

touch versions.yml
"""
}

process PREPARE_BED {
label "process_single"
label "array100"
tag "${meta.id}"
afterScript "rm -rf TMP"
conda   "${moduleDir}/../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta4),path(gtf)
    tuple val(meta ),path(gene_names)
output:
    tuple val(meta),path("*.simple.bed"),emit:bed
    tuple val(meta),path("*.exon.bed.gz"),path("*.exon.bed.gz.tbi"),path("*.exon.hdr"),emit:tabix
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:gtf.name
    def extend = 1000
"""
mkdir -p TMP
set -x

grep -v '^#' ${gene_names} | sort | uniq > TMP/genes.txt

jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP gtf2bed  \\
    --columns "gtf.feature,gene_name" -R "${fasta}"  "${gtf}" |\\
  awk -F '\t' '\$4=="exon" && \$5!="." && \$5!=""' |\\
  grep -f TMP/genes.txt -F |\\
  cut -f1,2,3,5 |\\
  sort --buffer-size=${task.memory.mega}M -t '\t' -k4,4 -T TMP  |\\
  uniq > TMP/exons.bed


join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,1.4' TMP/exons.bed TMP/genes.txt |\\
  sort -T TMP -t '\t' -k1,1 -k2,2n |\\
  uniq > ${prefix}.simple.bed


cat ${prefix}.simple.bed  | bgzip > ${prefix}.exon.bed.gz
tabix -p bed -f ${prefix}.exon.bed.gz

echo '##INFO=<ID=EXON,Number=.,Type=String,Description="EXON">' > ${prefix}.exon.hdr

bedtools slop -b ${extend} -g ${fai} -i ${prefix}.simple.bed |\\
  LC_ALL=C sort -t '\t' -k1,1 -k2,2n -T TMP > ${prefix}.extended.bed

touch versions.yml
"""
}


process ANNOT_SNV {
tag "${vcf.name}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
label "process_single"
label "array100"
afterScript "rm -rf TMP"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta4),path(bed)
    tuple val(meta ), path(vcf),path(vcfidx)
output:
    tuple val(meta ), path("*.vcf.gz"),path("*.vcf.gz.tbi"),optional:true,emit:vcf
    path("versions.yml"),emit:versions
script:
  def prefix = task.ext.prefix?:vcf.baseName+".annot"
  def gnomadAF = 0.001
  def gnomadPop = "AF_nfe"
  def soacn = "SO:0001818,SO:0001629"
  def jvm = "-Xmx${task.memory.giga}G  -Djava.io.tmpdir=TMP"
"""
hostname 1>&2
mkdir -p TMP

bcftools view --apply-filters '.,PASS' --regions-file "${bed}" "${vcf}" |\\
    bcftools norm  -f ${fasta} --multiallelics -any -O v |\\ 
    snpEff ${jvm} eff \\
        -dataDir "${params.snpeff_database_directory}" \\
        -nodownload -noNextProt -noMotif -noInteraction -noLog -noStats -chr chr -i vcf -o vcf "${params.snpeff_db}" |\\
    jvarkit ${jvm}  vcffilterso  --acn "${soacn}" |\\
    jvarkit ${jvm} vcfgnomad  --bufferSize 10000 \\
            --gnomad "${params.gnomad}" \\
            --fields "${gnomadPop}" \\
            --max-af "${gnomadAF}"  |\\
    bcftools view --apply-filters '.,PASS' > TMP/jeter.vcf


bcftools query -f '.' TMP/jeter.vcf > TMP/count.txt

if test -f TMP/count.txt
then
    bcftools query -l TMP/jeter.vcf > TMP/jeter.a
    bcftools query -l TMP/jeter.vcf | md5sum | cut -d ' ' -f 1 > TMP/jeter.b
    paste  TMP/jeter.a  TMP/jeter.b >  TMP/jeter.c

    bcftools reheader --samples TMP/jeter.c -T TMP/x  TMP/jeter.vcf  |\\
       bcftools view -O z -o "${prefix}.vcf.gz" 
    bcftools index -f -t "${prefix}.vcf.gz"
fi

touch versions.yml
"""
}



process ANNOT_SV {
tag "${vcf.name}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
label "process_single"
label "array100"
afterScript "rm -rf TMP"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta4),path(tabix),path(tabixidx),path(hdr)
    tuple val(meta5),path(gnomad_vcf),path(gnomad_vcf_tbi)
    tuple val(meta ), path(vcf),path(vcfidx)
output:
    tuple val(meta ), path("*.vcf.gz"),path("*.vcf.gz.tbi"),optional:true,emit:vcf
    path("versions.yml"),emit:versions
script:
  def prefix = task.ext.prefix?:vcf.baseName+".annot"
  def gnomadAF = 0.001
  def gnomadPop = "POPMAX_AF"
  def jvm = "-Xmx${task.memory.giga}G  -Djava.io.tmpdir=TMP"
"""
hostname 1>&2
mkdir -p TMP


bcftools view --apply-filters '.,PASS' -e 'SVTYPE=="BND" || SVTYPE=="INS" ' --regions-file "${tabix}" "${vcf}" |\\
    java -jar ${jvm} ${HOME}/jvarkit.jar vcfgnomadsv \\
            --gnomad "${gnomad_vcf}" \\
            --population "${gnomadPop}" \\
            --max-af "${gnomadAF}"  |\\
    bcftools annotate -a "${tabix}" -c 'CHROM,POS,END,EXON' -O u -h ${hdr} |\\
    bcftools view --apply-filters '.,PASS' > TMP/jeter.vcf


bcftools query -f '.' TMP/jeter.vcf > TMP/count.txt

if test -f TMP/count.txt
then
    bcftools query -l TMP/jeter.vcf > TMP/jeter.a
    bcftools query -l TMP/jeter.vcf | md5sum | cut -d ' ' -f 1 > TMP/jeter.b
    paste  TMP/jeter.a  TMP/jeter.b >  TMP/jeter.c

    bcftools reheader --samples TMP/jeter.c -T TMP/x  TMP/jeter.vcf  |\\
       bcftools view -O z -o "${prefix}.vcf.gz" 
    bcftools index -f -t "${prefix}.vcf.gz"
fi

touch versions.yml
"""
}



process MERGE_VARIANTS {
tag "${meta.id}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
label "process_single"
afterScript "rm -rf TMP"
input:
    tuple val(meta ), path("VCF/*")
output:
    tuple val(meta ), path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
script:
  def jvm = "-Xmx${task.memory.giga}G  -Djava.io.tmpdir=TMP"
  def prefix  = "${meta.id}.digest"
"""
mkdir -p TMP
find VCF/ -name "*.vcf.gz" | sort | uniq > TMP/jeter.list
jvarkit ${jvm} vcfconcat --merge  TMP/jeter.list > TMP/jeter.vcf
jvarkit ${jvm} vcfmulti2one TMP/jeter.vcf --discard_hom_ref  --discard_no_call  TMP/jeter.vcf > TMP/jeter2.vcf 
mv TMP/jeter2.vcf  TMP/jeter.vcf

bcftools sort --max-mem ${task.memory.giga}G -T TMP/sort  -O u TMP/jeter.vcf |\\
    bcftools view -O z -o "${prefix}.vcf.gz"
bcftools index -f -t "${prefix}.vcf.gz"
"""
}
