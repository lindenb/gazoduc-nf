
include {TRUVARI                                  } from '../../subworkflows/truvari/main.nf'
include {BCTOOLS_SETGT  as NO_CALL_TO_HOM_REF     } from '../../modules/bcftools/setgt/main.nf'
include {GATK_POSSIBLE_DENOVO  as   GATK_POSSIBLE_DENOVO_SV    } from '../../modules/gatk/possibledenovo/main.nf'
include {ANNOTATE_SV                              } from '../../subworkflows/annotation/sv/main.nf'


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
    TRUVARI(
	        meta,
            fasta,
            fai,
            dict,
            vcf
        )

    NO_CALL_TO_HOM_REF(TRUVARI.out.vcf)

    GATK_POSSIBLE_DENOVO_SV(fasta,fai,dict,pedigree,NO_CALL_TO_HOM_REF.out.vcf)
    FILTER_SV(fasta,fai,dict, GATK_POSSIBLE_DENOVO_SV.out.vcf)

    ANNOTATE_SV(
        meta,
        fasta,
        fai,
        dict,
        gtf,
        FILTER_SV.out.vcf
    )
        
emit:
    versions
}


process FILTER_SV {
tag "${meta.id?:""}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
    tuple val(meta),path(fasta)
    tuple val(meta),path(fai)
    tuple val(meta),path(dict)
    tuple val(meta),path(vcf),path(idx)
output:
    tuple val(meta),path("*.bcf"),path("*.csi"),emit:vcf
script:
	def prefix = task.ext.prefix?:meta.id+".filter"
"""
mkdir -p TMP
cat << EOF > TMP/jeter.code
String svType = variant.getAttribute("SVTYPE");
if(svType.equals("BND")) return false;

if(!(variant.hasAttribute("hiConfDeNovo") || variant.hasAttribute("loConfDeNovo))) return false;

int svlen ;
if(variant.hasAttribute("SVLEN")) {
    svlen = variant.getAttributeAsInt("SVLEN",0);
    }
else
    {
    svlen = variant.getEnd() - variant.getStart() +1;
    }
if(svlen<0) svlen = svlen*-1;

if(svlen < 500 ) return false;
if(svlen > 10000) return false;
return true;
EOF

bcftools view ${vcf} |\\
    jvarkit -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP vcffilterjdk \\
        -f  TMP/jeter.code |\\
    bcftools view --write-index -O b -o TMP/jeter.bcf

mv TMP/jeter.bcf ${prefix}.bcf
mv TMP/jeter.bcf.csi ${prefix}.bcf.csi
"""
}
