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
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

include {DOWNLOAD_GTF_OR_GFF3 as DOWNLOAD_GFF3    } from '../../modules/gtf/download/main.nf'
include {DOWNLOAD_GTF_OR_GFF3 as DOWNLOAD_GTF     } from '../../modules/gtf/download/main.nf'
include {SCATTER_TO_BED                           } from '../../subworkflows/gatk/scatterintervals2bed'
include {BCTOOLS_MENDELIAN2                       } from '../../modules/bcftools/mendelian2/main.nf'
include {HET_COMPOSITE  as APPPLY_HET_COMPOSITE   } from '../../modules/jvarkit/hetcomposite/main.nf'
include {BCFTOOLS_CONCAT as CONCAT1               } from '../../modules/bcftools/concat/main.nf'
include {WORKFLOW_SNV                             } from './sub.snv.nf'
include {WORKFLOW_SV                              } from './sub.sv.nf'
include {SOMALIER_BAMS                            } from '../../subworkflows/somalier/bams/main.nf'
include {MULTIQC                                  } from '../../modules/multiqc'
include {COMPILE_VERSIONS                         } from '../../modules/versions/main.nf'
include {runOnComplete; dumpParams                } from '../../modules/utils/functions.nf'
include { parseBoolean                            } from '../../modules/utils/functions.nf'
include { PREPARE_ONE_REFERENCE                   } from '../../subworkflows/samtools/prepare.one.ref'
include { VCF_INPUT                               } from '../../subworkflows/nf/vcf_input'
include { META_TO_BAMS                            } from '../../subworkflows/samtools/meta2bams2'
include { READ_SAMPLESHEET                        } from '../../subworkflows/nf/read_samplesheet'


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


boolean isStructuralVariantVCF(vcf) {
    if(vcf.name.endsWith(".sv.vcf.gz")) return true;//DRAGEN
    int found=0;
    try(InputStream in0= java.nio.file.Files.newInputStream(vcf)) {
        try(GZIPInputStream in=new GZIPInputStream(in0)) {
            try(BufferedReader br=new BufferedReader(new InputStreamReader(in))) {
                String line;
                while((line=br.readLine())!=null) {
                    if(!line.startsWith("#")) break;
                    if(line.equals("##source=Graphtyper")) return false;
                    if(line.equals("##source=DRAGEN_SV")) return true;
                    if(line.startsWith("##INFO=<ID=SVTYPE,")) found++;
                    if(line.startsWith("##INFO=<ID=SVLEN,")) found++;
                }
            }
        }
    }
    catch(Throwable err) {
	err.printStackTrace();
	throw err;
	}
    return found==2;
}

workflow {
        versions = Channel.empty()
        multiqc = Channel.empty()

        if(params.fasta==null) {
            log.error("--fasta missing")
            exit -1;
            }

         if(params.samplesheet==null) {
            log.error("--samplesheet missing")
            exit -1;
            }
        if(params.pedigree==null) {
            log.error("--pedigree missing")
            exit -1;
            }

        def metadata = [
            id: "trios",
            with_delly : params.with_delly,
            with_cardiopanel : params.with_cardiopanel , 
            with_alphamissense : params.with_alphamissense,
            with_no_avada : params.with_no_avada,
            with_clinvar : params.with_clinvar,
            with_snpeff : params.with_snpeff,
            with_vista : params.with_vista,
            with_vep : params.with_vep,
            with_utr_annotator : params.with_utr_annotator,
            with_vista : params.with_vista
            ];
        
        
        def pedigree = [[id:"pedigree"],file(params.pedigree)]

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

        /***************************************************
        *
        *  samplesheet contains VCFs
        *
        */
        VCF_INPUT(
            metadata.plus(
                path: params.samplesheet,
                require_index: true,
                arg_name : "samplesheet",
                required: true,
                unique : false
                )
            )
        versions = versions.mix(VCF_INPUT.out.versions)
        vcfs = VCF_INPUT.out.vcf


        bams = Channel.empty()
        triosbams_ch = Channel.empty()
        
        if(params.bams!=null) {
            READ_SAMPLESHEET(
                [arg_name:"bams"],
                params.bams
                )
             versions = versions.mix(READ_SAMPLESHEET.out.versions)

            META_TO_BAMS(
                metadata.plus([:]),
                PREPARE_ONE_REFERENCE.out.fasta,
                PREPARE_ONE_REFERENCE.out.fai,
                PREPARE_ONE_REFERENCE.out.dict,
                READ_SAMPLESHEET.out.samplesheet
                )
            versions = versions.mix(META_TO_BAMS.out.versions)


            bams = META_TO_BAMS.out.bams
            
            if(parseBoolean(params.with_somalier)) {
                SOMALIER_BAMS(
                    [id:"somalier"],
                    fasta,
                    fai,
                    dict,
                    META_TO_BAMS.out.map{meta,bam,bai,_fa,_fai,_dict->[meta,bam,bai]}, // sample,bam,bai
                    pedigree, // pedigree for somalier
                    Channel.of([[id:"no_sites"],[]])
                    )
                versions = versions.mix(SOMALIER_BAMS.out.versions)
                }

            triosbams_ch = Channel.fromPath(params.pedigree)
                .splitCsv(header:false,sep:'\t')
                .filter{!it[2].equals("0") && !it[3].equals("0")}
                .map{[it[1],it[2],it[3]]} // child,father,modther
                .combine(bams) // join child [ child,father,mother, metaC, bamC, baiC, fastaC, faiC, dictC ]
                .filter{it[0].equals(it[3].id)}
                .combine(bams) // join child [ child,father,mother, metaC, bamC, baiC, fastaC, faiC, dictC,  metaF, bamF, baiF, fastaF, faiF, dictF, ]
                .filter{it[1].equals(it[9].id) && it[6].equals(it[12])} //father and same fasta
                .combine(bams) // join child [ child,father,mother, metaC, bamC, baiC, fastaC, faiC, dictC,  metaF, bamF, baiF, fastaF, faiF, dictF, metaM, bamM, baiM, fastaM, faiM, dictM ]
                .filter{it[2].equals(it[15].id) && it[6].equals(it[18])} //father and same fasta
                .map{[
                    it[0], //child
                    it[1], //father
                    it[2], //mother
                    it[6], //fasta
                    it[7], //fai
                    it[8], //dict
                    it[4], it[5], // bamC, baiC
                    it[10],it[11], // bamP, baiP
                    it[16],it[17] // bamM, baiM
                ]}
        }   

      

        vcfs = vcfs.branch{
            sv : isStructuralVariantVCF(it[1])
            snv: true
         }

        DOWNLOAD_GFF3(PREPARE_ONE_REFERENCE.out.dict)
        versions = versions.mix(DOWNLOAD_GFF3.out.versions)

        DOWNLOAD_GTF(PREPARE_ONE_REFERENCE.out.dict)
        versions = versions.mix(DOWNLOAD_GTF.out.versions)

        if(params.bed==null) {
            bed = PREPARE_ONE_REFERENCE.out.scatter_bed
         } else {
             bed = [metadata, file(params.bed)]
         }



        WORKFLOW_SNV(
            metadata.plus([id:"triosnv"]),
            PREPARE_ONE_REFERENCE.out.fasta,
            PREPARE_ONE_REFERENCE.out.fai,
            PREPARE_ONE_REFERENCE.out.dict,
            bed,
            pedigree,
            [metadata,file(params.gnomad),file(params.gnomad+".tbi")],
            DOWNLOAD_GFF3.out.gtf,
            DOWNLOAD_GTF.out.gtf,
            triosbams_ch,
            vcfs.snv
            )
        
        
        versions = versions.mix(WORKFLOW_SNV.out.versions)
        multiqc =  multiqc.mix(WORKFLOW_SNV.out.multiqc)
        
        /*
        WORKFLOW_SV(
            metadata,
            PREPARE_ONE_REFERENCE.out.fasta,
            PREPARE_ONE_REFERENCE.out.fai,
            PREPARE_ONE_REFERENCE.out.dict,
            bed,
            pedigree,
            DOWNLOAD_GFF3.out.gtf,
            DOWNLOAD_GTF.out.gtf,
            triosbams_ch,
            vcfs.sv
            )
            
        versions = versions.mix(WORKFLOW_SV.out.versions)
        */
        
/*
        COMPILE_VERSIONS(versions.collect())
        multiqc = multiqc.mix(COMPILE_VERSIONS.out.multiqc)

        MULTIQC(multiqc.collect().map{[[id:"trios"],it]})
  */     
}

runOnComplete(workflow)


process README {
tag "${meta.id?:""}"
input:
    val(meta)
output:
    path("*.md"),emit:readme
script:
"""
cat << _EOF_ > README.md

The directory SNV_DENOVO contains

- a VCF file *.denovo.vcf.gz  containing all the "protein-altering variants" variants containing the de-novo variants (including low and high quality). The VCF was annotated with SNPEFF, VEP, BCFTOOLS_CSQ, CADD, CLINVAR, ALPHAMISSENSE, REVEL, the CARDIOPANEL list,  REPEATMASKER, GNOMAD, etc... The "INFO/CONTROLS_HAVING_ALT" is a count of parents from other trios carrying a ALT allele and could be used to remove false positives. The de-novo variants have been extracted using "GATK possibledeNovo" (INFO/hiConfDeNovo and loConfDeNovo).

- *.denovo.table.txt is a "vertical" user-friendly, view of the variant.
- *.denovo.genes.tsv is a summary by gene

The  directory 'results/IGV' contains the IGV profiles of the "hiConfidence de novo variant". Some of them are nevertheless obviously **FALSE POSITIVES**

The directory results/HETCOMPOSITE contains the  protein-altering variants, ( fitered on gnomad AF_NFE=1%) that could be used as pairs of het-composites variants.

```
results/
|-- HETCOMPOSITE
|   |-- *.concat.hetcomposite.vcf.gz
|   |-- *.concat.hetcomposite.vcf.gz.tbi
|   |-- *.concat.hetcomposite.genes.report.bed
|   `-- *.concat.hetcomposite.variants.report.txt
|-- IGV
|   |-- ****-C.html
|   `-- index.html
`-- SNV_DENOVO
    |-- *.denovo.genes.tsv
    |-- *.denovo.table.txt
    |-- *.denovo.vcf.gz
    `-- *.denovo.vcf.gz.tbi
```
_EOF_
"""
}
