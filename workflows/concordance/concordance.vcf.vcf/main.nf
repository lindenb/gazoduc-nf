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
include { runOnComplete                 } from '../../../modules/utils/functions.nf'
include { PREPARE_ONE_REFERENCE         } from '../../../subworkflows/samtools/prepare.one.ref'
include { MULTIQC                       } from '../../../modules/multiqc'
include { COMPILE_VERSIONS              } from '../../../modules/versions/main.nf'
include { VCF_INPUT as VCF_INPUT_TRUTH  } from '../../../subworkflows/nf/vcf_input'
include { VCF_INPUT as VCF_INPUT_CALL   } from '../../../subworkflows/nf/vcf_input'
include { BED_CLUSTER                   } from '../../../modules/jvarkit/bedcluster'
include { BCFTOOLS_MERGE                } from '../../../modules/bcftools/merge'
include { BCFTOOLS_CONCAT               } from '../../../modules/bcftools/concat3'
include { BCFTOOLS_INDEX                } from '../../../modules/bcftools/index'
include { BCFTOOLS_QUERY                } from '../../../modules/bcftools/query'
include { BCFTOOLS_QUERY_SAMPLES        } from '../../../modules/bcftools/query.samples'
include { JVARKIT_VCFGNOMAD             } from '../../../modules/jvarkit/vcfgnomad'
include { GHOSTSCRIPT_MERGE             } from '../../../modules/gs/merge'

workflow {
	
	if(params.fasta==null) {
		log.warn("undefined --fasta");
		java.lang.System.exit(-1);
		}
	if(params.truth_vcf==null) {
		log.warn("undefined --truth_vcf");
		java.lang.System.exit(-1);
		}
    if(params.call_vcf==null) {
		log.warn("undefined --call_vcf");
		java.lang.System.exit(-1);
		}
    
   versions = Channel.empty()
   multiqc = Channel.empty()
   metadata = [
		id:"concordance"
		]

    PREPARE_ONE_REFERENCE(
   		metadata.plus(skip_scatter:true),
   		Channel.of(params.fasta).map{f->file(f)}.map{f->[[id:f.baseName],file(f)]}.first()
   		)
    versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)


    /* if it's an exome , group the small genome together in BED */
    input_bed_ch =  [[id:"bed"],file(params.bed)]
    BED_CLUSTER(
        PREPARE_ONE_REFERENCE.out.dict,
        input_bed_ch
        )
    versions = versions.mix(BED_CLUSTER.out.versions)
    bed_ch = BED_CLUSTER.out.bed
        .map{_meta,bed->bed}
        .map{it instanceof List?it:[it]}
        .flatMap()
        .map{bed->[[id:bed.baseName],bed]}
    
   VCF_INPUT_TRUTH(metadata.plus([
            path: params.truth_vcf,
            arg_name: "truth_vcf",
            require_index : true,
            required: true,
            unique : false
            ]))
    versions = versions.mix(VCF_INPUT_TRUTH.out.versions)


    all_truth_bed_ch = VCF_INPUT_TRUTH.out.vcf
        .map{_meta,vcf,tbi->[vcf,tbi]}
        .flatMap()
        .collect()
        .map{files->[files.sort()]}
        .combine(bed_ch)
        .map{files,meta,bed->[meta.plus(type:"truth_vcf"),bed,files]}
        .multiMap{meta,bed,files->
            bed: [meta,bed]
            vcf: [meta, files]
        }

    BCFTOOLS_MERGE(
        all_truth_bed_ch.bed,
        all_truth_bed_ch.vcf
        )
    versions = versions.mix(BCFTOOLS_MERGE.out.versions)

    gnomad_ch = [
        [id:"gnomad"],
        file(params.gnomad),
        file(params.gnomad+".tbi")
        ]

    JVARKIT_VCFGNOMAD(
        gnomad_ch,
        BCFTOOLS_MERGE.out.vcf.map{meta,vcf,_tbi->[meta,vcf]}
        )
    versions = versions.mix(JVARKIT_VCFGNOMAD.out.versions)

    BCFTOOLS_INDEX(JVARKIT_VCFGNOMAD.out.vcf)
    versions = versions.mix(BCFTOOLS_INDEX.out.versions)

    BCFTOOLS_CONCAT(
        BCFTOOLS_INDEX.out.vcf  
        .map{_meta,vcf,tbi->[vcf,tbi]}
        .flatMap()
        .collect()
        .map{files->[[id:"truth"],files.sort()]}
        )
    versions = versions.mix(BCFTOOLS_CONCAT.out.versions)


  /***************************************************
    *
    * VCF input
    *
    */
    VCF_INPUT_CALL(metadata.plus([
            path: params.call_vcf,
            arg_name: "call_vcf",
            require_index : true,
            required: true,
            unique : false
            ]))
    versions = versions.mix(VCF_INPUT_CALL.out.versions)
    
    CONCORDANCE(
        PREPARE_ONE_REFERENCE.out.fasta,
        PREPARE_ONE_REFERENCE.out.fai,
        PREPARE_ONE_REFERENCE.out.dict,
        BCFTOOLS_CONCAT.out.vcf,
        VCF_INPUT_CALL.out.vcf
        )
    versions = versions.mix(CONCORDANCE.out.versions)

    PLOT_CONCORDANCE(CONCORDANCE.out.concordance)
    versions = versions.mix(PLOT_CONCORDANCE.out.versions)

    MERGE_CONCORDANCES( 
		CONCORDANCE.out.concordance
			.map{_meta,f->f}
			.collect()
			.map{files->[metadata,files.sort()]}
		)
	versions = versions.mix(MERGE_CONCORDANCES.out.versions)


	GHOSTSCRIPT_MERGE(
		PLOT_CONCORDANCE.out.pdf
			.map{_meta,f->f}
			.collect()
			.map{files->[ [id:"plot"],files.sort()]}
		)
   	versions = versions.mix(GHOSTSCRIPT_MERGE.out.versions)


    COMPILE_VERSIONS(versions.collect())
    multiqc = multiqc.mix(COMPILE_VERSIONS.out.multiqc)


    }

runOnComplete(workflow)

process CONCORDANCE {
tag "${meta.id}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
    tuple val(meta6),path(truth_vcf),path(truth_vcf_tbi)
    tuple val(meta ),path(call_vcf),path(tbi)
output:
    tuple val(meta),path("*.txt.gz"),emit:concordance
    path("versions.yml"),emit:versions
script:
  def prefix = task.ext.prefix?:"${meta.id}"
  def jvm = task.ext.jvm?:" -Djava.io.tmpdir=TMP"
"""
mkdir -p TMP/isec1  TMP/isec2
set -x
bcftools isec -c none ${call_vcf} "${truth_vcf}" -O z  -p TMP/isec1 -n =2 -w 1
find TMP/isec1 -type f 1>&2
cat TMP/isec1/sites.txt 1>&2

cat << 'EOF' > TMP/jeter.code
if(!variant.getContig().matches("(chr)?[0-9]+")) return false;
if(variant.isFiltered()) return false;

for(Genotype g : variant.getGenotypes()) {
    if(g.isFiltered()) return false;
    if(!g.isHet()) return false;
    if(!g.hasAD()) return false;
    final int[] ad = g.getAD();
    if(ad.length!=2) return false;
    if(ad[0]< 2 || ad[1]<2) return false;
    return true;
	}
return false;
EOF

bcftools view  --apply-filters 'PASS,.' -m2 -M2 --types snps -Ov  TMP/isec1/0000.vcf.gz |\\
   jvarkit ${jvm} vcffilterjdk -f TMP/jeter.code  |\\
   bcftools view -O z -o TMP/jeter2.vcf.gz

bcftools index -t -f  TMP/jeter2.vcf.gz

bcftools isec --write-index -c none TMP/jeter2.vcf.gz "${truth_vcf}" -O b  -p TMP/isec2 -n =2 -w 2
find TMP/isec2 -type f  1>&2
cat TMP/isec2/sites.txt 1>&2

bcftools gtcheck -u GT -g TMP/jeter2.vcf.gz  TMP/isec2/0001.bcf |\
	grep  '^DCv2' | cut -f2- |\\
	awk '{printf("%s\t%f\\n",\$0,(int(\$5)==0?-1:(int(\$6)/(1.0*int(\$5)))));}'  |\
	LC_ALL=C sort -T TMP -t '\t' -k7,7gr |\\
	awk -F '\t' '{printf("%s\t%s\\n",\$0,(NR==1 && \$7*1.0 > 0.0?"best":"."));}' |gzip > ${prefix}.concordance.txt.gz

touch versions.yml
"""
}



process PLOT_CONCORDANCE {
label "process_single"
tag "${meta.id}"
conda "${moduleDir}/../../../conda/bioinfo.02.yml"
afterScript "rm -rf TMP"
input:
        tuple val(meta),path(metrics)
		
output:
        tuple val(meta),path("*.pdf"),emit:pdf
		path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:"${meta.id}.concordances"
	def treshold = (task.ext.treshold?:0.7) as double
"""
mkdir -p TMP
gunzip -c '${metrics}' > TMP/jeter.tsv
cat << '__EOF__' | R --vanilla --no-save
table <- read.table("TMP/jeter.tsv",sep="\t",header=FALSE,stringsAsFactors=FALSE);
colnames(table) <-c("S1","S2","x1","x2","x3","GT_ALT","CONCORDANCE","x5")

maxc = max(table\$CONCORDANCE)

pdf(paste0(((1.0-maxc)*100.0),"${prefix}.pdf"));
plot(
        x = table\$GT_ALT,
        y = table\$CONCORDANCE,
        xlab="Number_of_matching_genotypes",
        ylab="Concordance",
        ylim=c(0,1.0),
        main="Concordances ${meta.id}",
        sub=""
        );


best <- table[table\$CONCORDANCE == maxc | table\$CONCORDANCE > ${treshold},]
text(x = best\$GT_ALT, y= best\$CONCORDANCE, labels = best\$S1 , pos = 2 , cex = 0.5 , col= "green") 

abline(h=${treshold},col="blue")
dev.off();
__EOF__

touch versions.yml
"""
stub:
	def prefix = task.ext.prefix?:"${meta.id}.concordances"
"""
touch versions.yml ${prefix}.pdf
"""
}




process MERGE_CONCORDANCES {
label "process_single"
tag "${meta.id}"
afterScript "rm -rf TMP"
input:
	tuple val(meta),path("DATA/??.txt")
output:
	tuple val(meta),path("*.tsv.gz"),emit:table
	path("versions.yml"),emit:versions
script:
	def prefix= task.ext.prefix?:"${meta.id}.concordance"
"""
mkdir -p TMP
echo -e 'Query_Sample\tGenotyped_Sample\tDiscordance\tAverage_-log_P(HWE)\tNumber_of_sites_compared\tNumber_of_matching_genotypes\tfraction\tbest' |gzip > ${prefix}.tsv.gz
find DATA -name "*.txt" -exec gunzip -c '{}' ';'  | LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k7,7gr -k6,6nr | gzip --best >> ${prefix}.tsv.gz
touch versions.yml
"""
stub:
	def prefix= task.ext.prefix?:"${meta.id}.concordance"
"""
touch versions.yml ${prefix}.tsv.gz
"""
}
