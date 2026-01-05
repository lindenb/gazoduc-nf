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

include { BED_STATS                  } from '../../modules/bedstats'
include { dumpParams;runOnComplete   } from '../../modules/utils/functions.nf'
include { SCATTER_TO_BED             } from '../../subworkflows/gatk/scatterintervals2bed'
include { MANTA_MULTI                } from '../../modules/manta/multi'
include { MANTA_CONVERT_INVERSION    } from '../../modules/manta/convert.inversion'
include { TRUVARI_COLLAPSE           } from '../../modules/truvari/collapse'
include { ANNOTATE_SV                } from '../../subworkflows/annotation/sv'
include { JVARKIT_VCF_SET_DICTIONARY } from '../../modules/jvarkit/vcfsetdict'
include { MULTIQC                    } from '../../modules/multiqc'
include { COMPILE_VERSIONS           } from '../../modules/versions/main.nf'


if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}


Map assertKeyExists(final Map hash,final String key) {
    if(!hash.containsKey(key)) throw new IllegalArgumentException("no key ${key}'in ${hash}");
    return hash;
}

Map assertKeyExistsAndNotEmpty(final Map hash,final String key) {
    assertKeyExists(hash,key);
    def value = hash.get(key);
    if(value.isEmpty()) throw new IllegalArgumentException("empty ${key}'in ${hash}");
    return hash;
}

Map assertKeyMatchRegex(final Map hash,final String key,final String regex) {
    assertKeyExists(hash,key);
    def value = hash.get(key);
    if(!value.matches(regex)) throw new IllegalArgumentException(" ${key}'in ${hash} doesn't match regex '${regex}'.");
    return hash;
}

workflow {

    def ref_hash = [
        id: file(params.fasta).simpleName,
        ucsc_name: (params.ucsc_name?:"undefined"),
        ensembl_name: (params.ensembl_name?:"undefined"),
        ]
    def fasta      = [ref_hash, file(params.fasta) ]
    def fai        = [ref_hash, file(params.fai) ]
    def dict       = [ref_hash, file(params.dict) ]
    def gtf        = [ref_hash, file(params.gtf), file(params.gtf+".tbi") ]
    def pedigree   = [[id:"nopedigree"], [] ]

    if(params.pedigree!=null) {
         pedigree   = [[id:file(params.pedigree).baseName], file(params.pedigree) ]
    }

    versions = Channel.empty()
    multiqc = Channel.empty()

    ch1 = Channel.fromPath(params.samplesheet)
        .splitCsv(header:true,sep:',')
        .map{assertKeyMatchRegex(it,"sample","^[A-Za-z_0-9\\.\\-]+\$")}
        .map{assertKeyMatchRegex(it,"bam","^\\S+\\.(bam|cram)\$")}
        .map{
            if(it.containsKey("batch")) return it;
            return it.plus(batch:"batch_"+it.sample);
        }
        .map{assertKeyMatchRegex(it,"batch","^[A-Za-z_0-9\\.\\-]+\$")}
        .map{
            if(it.containsKey("bai")) return it;
            if(it.bam.endsWith(".cram")) return it.plus(bai : it.bam+".crai");
            return it.plus(bai:it.bam+".bai");
        }
        .map{
            if(it.containsKey("fasta")) return it;
            return it.plus(fasta:params.fasta);
        }
        .map{
            if(it.containsKey("fai")) return it;
            return it.plus(fai:it.fasta+".fai");
        }
        .map{
            if(it.containsKey("dict")) return it;
            return it.plus(dict: it.fasta.replaceAll("\\.(fasta|fa|fna)\$",".dict"));
        }
        .map{assertKeyMatchRegex(it,"bai","^\\S+\\.(bai|crai)\$")}
        .map{[ it.batch, it.sample,file(it.bam),file(it.bai),file(it.fasta),file(it.fai),file(it.dict)]}
        .groupTuple()
        .map{
            java.util.HashSet<String> set = new  java.util.HashSet<String>();
            set.addAll(it[4]);
            if(set.size()!=1) throw new IllegalArgumentException("multiple fasta file for ${it} for the same batch.")
            return it;
        }
        .multiMap {
            fasta: [ref_hash,it[4][0]]
            fai: [ref_hash,it[5][0]]
            dict: [ref_hash,it[6][0]]
            bams: [ref_hash.plus(id:it[0]), it[2].plus(it[3]) ]
        }
   
                
    if(params.bed==null) {
        SCATTER_TO_BED(ref_hash,fasta,fai,dict)
        versions = versions.mix(SCATTER_TO_BED.out.versions)
        bed = SCATTER_TO_BED.out.bed
        }
	else {
        bed =  [ref_hash, file(params.bed) ]
        }


    BED_STATS(fai,bed)
    versions = versions.mix(BED_STATS.out.versions)
    multiqc = multiqc.mix(BED_STATS.out.multiqc)

    MANTA_MULTI(
        ch1.fasta,
        ch1.fai,
        ch1.dict,
        bed,
        ch1.bams
        )
    versions = versions.mix(MANTA_MULTI.out.versions.first())
    
    JVARKIT_VCF_SET_DICTIONARY(
        dict,
        MANTA_MULTI.out.diploidSV
        )
    versions = versions.mix(JVARKIT_VCF_SET_DICTIONARY.out.versions.first())


    MERGE_CANDIDATE_SV(
        dict,
        MANTA_MULTI.out.candidateSV .map{[it[1],it[2]]}
            .collect()
            .map{[[id:"candidateSV"],it.flatten()]}
        )

    MANTA_CONVERT_INVERSION(
        fasta,
        fai,
        dict,
        JVARKIT_VCF_SET_DICTIONARY.out.vcf
        )
    versions = versions.mix(MANTA_CONVERT_INVERSION.out.versions.first())

    if(params.pedigree!=null) {
        MANTA_DENOVO(
            fasta,
            fai,
            dict,
            pedigree,
            MANTA_CONVERT_INVERSION.out.vcf
            )
        versions = versions.mix(MANTA_DENOVO.out.versions)
        }



    TRUVARI_COLLAPSE(
        fasta,
        fai,
        dict,
        MANTA_CONVERT_INVERSION.out.vcf
            .map{[it[1],it[2]]}
            .collect()
            .map{[[id:"truvari"],it.flatten()]}
        )
    versions = versions.mix(TRUVARI_COLLAPSE.out.versions.first())
	
    ANNOTATE_SV(
        ref_hash,
        fasta,
        fai,
        dict,
        gtf,
        TRUVARI_COLLAPSE.out.vcf
        )
    versions = versions.mix(ANNOTATE_SV.out.versions.first())


    COMPILE_VERSIONS(versions.collect())
    multiqc = multiqc.mix(COMPILE_VERSIONS.out.multiqc)

    MULTIQC([[id:"no_mqc_config"],[]],
        multiqc.collect().map{[[id:"manta"],it]})
    }


runOnComplete(workflow)

process MERGE_CANDIDATE_SV {
    label "process_single"
    tag "${meta.id?:""}"
    afterScript "rm -rf TMP"
    conda "${moduleDir}/../../conda/bioinfo.01.yml"
    input:
		tuple val(meta3),path(dict)
        tuple val(meta),path("VCFS/*")
    output:
        tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
    	path("versions.yml"),emit:versions
    script:
        def prefix = task.ext.prefix?:meta.id + ".candidateSV"
    """
    mkdir -p TMP
    find VCFS/ -name "*.vcf.gz" > TMP/jeter.list

    bcftools concat --allow-overlaps -O v --drop-genotypes --rm-dups exact  --file-list  TMP/jeter.list |\\
            jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP vcfsetdict \\
                    -n SKIP \\
                    --reference '${dict}' |\\
            bcftools annotate -i 'SVTYPE!="BND"' --set-id "%CHROM:%POS:%END:%SVTYPE" -O u -o TMP/jeter.bcf

    bcftools view --header-only TMP/jeter.bcf > TMP/jeter2.vcf

    bcftools view --no-header TMP/jeter.bcf |\\
        sort -T TMP -t '\t' -k3,3 --unique |\\
        sort -T TMP -t '\t' -k1,1 -k2,2n >> TMP/jeter2.vcf


    bcftools view -O z -o TMP/jeter.vcf.gz TMP/jeter2.vcf
    bcftools index -f -t  TMP/jeter.vcf.gz


    mv  TMP/jeter.vcf.gz "${prefix}.vcf.gz"
    mv  TMP/jeter.vcf.gz.tbi "${prefix}.vcf.gz.tbi"

    cat << EOF > versions.yml
    "${task.process}":
        jvarkit: todo
    EOF
    """
    }


process MANTA_DENOVO {
    label "process_single"
    tag "${meta.id?:""} ${vcf.name}"
    afterScript "rm -rf TMP"
    conda "${moduleDir}/../../conda/bioinfo.01.yml"
    input:
		tuple val(meta1),path(fasta)
		tuple val(meta2),path(fai)
		tuple val(meta3),path(dict)
        tuple val(meta4),path(pedigree)
        tuple val(meta),path(vcf),path(vcfidx)
    output:
        tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
    	path("versions.yml"),emit:versions
  
    script:
        def prefix = task.ext.prefix?:meta.id + ".denovo"
        def gq = task.ext.gq?:50
	"""
	hostname 1>&2
	mkdir -p TMP

cat << __EOF__ > TMP/jeter.code

private boolean acceptControl(final VariantContext vc,String sm) {
        final Genotype g = vc.getGenotype(sm);
        if(g!=null && g.hasAltAllele()) return false;
        return true;
        }

private boolean acceptTrio(final VariantContext vc,String cm,String fm,String mm) {
    final Genotype c = vc.getGenotype(cm);
    if(c==null) return false;
    if(c.isFiltered() && !c.getFilters().equals("PASS")) return false;
    if(c.isHomVar()) return false;
    if(c.hasGQ() && c.getGQ() < ${gq} ) return false;
    final Genotype m = vc.getGenotype(mm);
    final Genotype f = vc.getGenotype(fm);
    // child must be HET
    if(!c.hasAltAllele()) return false;
    // keep de novo
    if(f.hasAltAllele() || m.hasAltAllele()) return false;
    if(!f.getFilters().matches("(HomRef|PASS)")) return false;
    if(f.hasGQ() && f.getGQ() < ${gq} ) return false;
    if(!m.getFilters().matches("(HomRef|PASS)")) return false;
    if(m.hasGQ() && m.getGQ() < ${gq} ) return false;
    return true;
    }

public Object apply(final VariantContext variant) {
if(variant.isFiltered()) return false;
final String svType = variant.getAttributeAsString("SVTYPE","");
if(svType.equals("BND")) return false;
__EOF__



    # convert pedigree if no 6th column
    awk '{S=\$6 ; if(NF==5 || S=="") { if(\$3!="0" && \$4!="0") {S="case";} else {S="control"} }  printf("%s\t%s\t%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,\$4,\$5,S);}' ${pedigree} > TMP/pedigree.tsv
    

    ## all other samples are controls
comm -13 \\
        <(cut -f 2 TMP/pedigree.tsv  | sort | uniq) \\
        <(bcftools query -l '${vcf}'| sort | uniq) |\\
        awk '{printf("if(!acceptControl(variant,\\"%s\\")) return false;\\n",\$1);}' >> TMP/jeter.code

awk -F '\t' '(\$6=="control" || \$6=="unaffected") {printf("if(!acceptControl(variant,\\"%s\\")) return false;\\n",\$2);}'  TMP/pedigree.tsv >> TMP/jeter.code

awk -F '\t' '(\$6=="case" || \$6=="affected") {printf("if(acceptTrio(variant,\\"%s\\",\\"%s\\",\\"%s\\")) return true;\\n",\$2,\$3,\$4);}' TMP/pedigree.tsv  >> TMP/jeter.code

cat << EOF >> TMP/jeter.code
return false;
}
EOF

    bcftools view "${vcf}"  |\\
            jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP vcffilterjdk --body -f TMP/jeter.code |\\
            bcftools view -o TMP/jeter.vcf.gz -O z

    bcftools index --threads ${task.cpus} -f -t TMP/jeter.vcf.gz


    mv  TMP/jeter.vcf.gz "${prefix}.vcf.gz"
    mv  TMP/jeter.vcf.gz.tbi "${prefix}.vcf.gz.tbi"


cat << EOF > versions.yml
"${task.process}":
	jvarkit: todo
EOF
"""
}


