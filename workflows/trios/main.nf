import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

include {CLINVAR                              } from '../../subworkflows/annotation/clinvar/main.nf'
include {ALPHAMISSENSE                        } from '../../subworkflows/annotation/alphamissense/main.nf'
include {VEP                                  } from '../../subworkflows/annotation/vep/main.nf'
include {TRIOS  as TRIO_SNV                   } from '../../subworkflows/trios/main.nf'
include {TRUVARI                              } from '../../subworkflows/truvari/main.nf'
include {SNPEFF                               } from '../../subworkflows/snpeff/main.nf'
include {DOWNLOAD_GENCODE as DOWNLOAD_GFF3    } from '../../modules/gtf/download/main.nf'
include {BCFTOOLS_BCSQ                        } from '../../modules/bcftools/bcsq/main.nf'
include {ANNOTATE_SV                          } from '../../subworkflows/annotation/sv/main.nf'
include {BCTOOLS_SETGT  as NO_CALL_TO_HOM_REF } from '../../modules/bcftools/setgt/main.nf'
include {BCTOOLS_MENDELIAN2                   } from '../../modules/bcftools/mendelian2/main.nf'


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
        def ref_hash = [
            id: file(params.fasta).simpleName
            ]
        def fasta  = [ref_hash, file(params.fasta) ]
        def fai    = [ref_hash, file(params.fai) ]
        def dict   = [ref_hash, file(params.dict) ]
        def NO_BED = [ref_hash, [] ]
        def pedigree = [[id:"pedigree"],file(params.pedigree)]

        vcfs = Channel.fromPath(params.samplesheet)
            .splitCsv(header:true,sep:',')
            .map{assertKeyMatchRegex(it,"vcf","^\\S+\\.(vcf\\.gz|bcf)\$")}
            .map{
                if(it.containsKey("idx")) return it;
                if(it.vcf.endsWith(".bcf")) return it.plus(idx : it.vcf+".csi");
                return it.plus(idx:it.vcf+".tbi");
            }
            .map{
                 if(it.containsKey("id")) return it;
                return it.plus(id:file(it.vcf).name.replaceAll("\\.(vcf\\.gz|bcf)\$",""));
                }
            .map{assertKeyMatchRegex(it,"idx","^\\S+\\.(tbi|csi)\$")}
            .map{[[id:it.id],file(it.vcf),file(it.idx)]}


        vcfs = vcfs.branch{
            sv : isStructuralVariantVCF(it[1])
            snv: true
         }


       DOWNLOAD_GFF3(fasta,fai,dict)


        


         STRUCTURAL_VARIANTS(
            [id:"sv"],
            fasta,
            fai,
            dict,
            DOWNLOAD_GFF3.out.output,
            pedigree,
            vcfs.sv.take(0)
         )


        vcf2bed = VCF_TO_CONTIGS(vcfs.snv) 
        vcf_snv = vcf2bed.output
            .splitText(elem:1)
            .map{[it[0],it[1].trim(),it[2],it[3]]}
            .map{[it[0].plus(contig:it[1]),it[2],it[3]]}

        
        
        SUBSET_VCF(fasta,fai, vcf_snv)
        vcf_snv = SUBSET_VCF.out.vcf

        TRIO_SNV(
            [id:"trio"],
            fasta,
            fai,
            dict,
            pedigree,
            vcf_snv
            )
        vcf_snv = TRIO_SNV.out.vcf

        
        VEP([id:"vep"],fasta,fai,dict,vcf_snv)
        vcf_snv = SUBSET_VCF.out.vcf

        SNPEFF([id:"snpeff"],fasta,fai,dict,vcf_snv)
        vcf_snv = SNPEFF.out.vcf

        BCFTOOLS_BCSQ(fasta,fai,DOWNLOAD_GFF3.out.output,vcf_snv )
        vcf_snv = BCFTOOLS_BCSQ.out.vcf

        FILTER_PREDICTIONS(vcf_snv)
        vcf_snv = FILTER_PREDICTIONS.out.vcf //we can use this later for het composite

        FILTER_MENDELIAN(vcf_snv)
        vcf_snv = SUBSET_VCF.out.vcf

        vcf_snv  = Channel.empty()

        CLINVAR([id:"clinvar"],fasta,fai,dict,NO_BED,vcf_snv)
        vcf_snv = CLINVAR.out.vcf
 
        ALPHAMISSENSE([id:"alphamissense"],fasta,fai,dict,NO_BED,vcf_snv)
        vcf_snv = ALPHAMISSENSE.out.vcf
        
        
     
	
        MERGE_AND_FILTER(
            vcf_snv
                .map{[[id:"snv"],it[1],it[2]]}
                .groupTuple()
                .map{[it[0],it[1].flatten()]}
        )

}
process VCF_TO_CONTIGS{
    tag "${meta.id}"
    label "process_single"
    conda "${moduleDir}/../../conda/bioinfo.01.yml"
    input:
        tuple val(meta),path(vcf),path(idx)
    output:
        tuple val(meta),path("chroms.txt"),path(vcf),path(idx),emit:output
    script:
        """
        bcftools index -s "${vcf}" | cut -f1 > chroms.txt
        """
}

process SUBSET_VCF {
    tag "${meta.id} ${meta.contig?:""}"
    label "process_single"
    conda "${moduleDir}/../../conda/bioinfo.01.yml"
    input:
        tuple val(meta1),path(fasta)
        tuple val(meta2),path(fai)
        tuple val(meta ),path(vcf),path(idx)
    output:
        tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
    script:
        def args = meta.contig?"--regions \"${meta.contig}\"":""
        def prefix = task.ext.prefix?:meta.id+(meta.contig?"."+meta.contig:"")
        """
        bcftools norm --multiallelics -both --fasta-ref "${fasta}" -O u  ${args}  "${vcf}" |\\
        bcftools view -e 'ALT=="*"' -Ou |\\
        bcftools annotate --set-id '%VKX'  -O z -o "${prefix}.vcf.gz"
        bcftools index --threads "${task.cpus}" -t -f "${prefix}.vcf.gz"
        """
}



process FILTER_MENDELIAN {
    tag "${meta.id?:vcf.name}"
    label "process_single"
    conda "${moduleDir}/../../conda/bioinfo.01.yml"
    input:
        tuple val(meta ),path(vcf),path(idx)
    output:
        tuple val(meta),path("*.bcf"),path("*.csi"),emit:vcf
    script:
        def prefix = task.ext.prefix?:vcf.baseName+".merr"
    """
    bcftools view --threads "${task.cpus}" -i 'INFO/loConfDeNovo!="" || INFO/hiConfDeNovo!="" || INFO/MERR>0' -O b -o ${prefix}.bcf '${vcf}'
    bcftools index --threads "${task.cpus}"  -f "${prefix}.bcf"
    """
}

process FILTER_PREDICTIONS {
    tag "${meta.id}"
    label "process_single"
    conda "${moduleDir}/../../conda/bioinfo.01.yml"
    afterScript "rm -rf TMP"
    input:
        tuple val(meta ),path(vcf),path(idx)
    output:
        tuple val(meta),path("*.bcf"),path("*.csi"),emit:vcf
    script:
        def args = meta.contig?"--regions \"${meta.contig}\"":""
        def prefix = task.ext.prefix?:vcf.baseName+".filterso"
        def acn = task.ext.acn?:"SO:0001818,SO:0001629"
        """
        mkdir -p TMP
        bcftools view "${vcf}" |\\
            jvarkit vcffilterso  -A '${acn}' --filterout BAD_SO |\\
            bcftools view  -O b -o "${prefix}.vcf.gz"
        bcftools index --threads "${task.cpus}" -t -f "${prefix}.vcf.gz"
        """
}

process MERGE_AND_FILTER {
    tag "${meta.id} ${meta.contig?:""}"
    label "process_single"
    conda "${moduleDir}/../../conda/bioinfo.01.yml"
    afterScript "rm -rf TMP"
    input:
        tuple val(meta),path("VCFS/*")
    script:
        def prefix = task.ext.prefix?:"snv"
    """
    mkdir -p TMP
cat << EOF > TMP/jeter.code
return false;
EOF

    find VCFS/ -name "*.vcf.gz" -o -name "*.bcf" > TMP/jeter.list
    
    bcftools concat --allow-overlaps --file-list TMP/jeter.list -O v |\\
        java -jar ~/jvarkit.jar vcffilterjdk -f TMP/jeter.code |\\
        bcftools view  -O z -o TMP/jeter.vcf.gz

    bcftools index -t -f TMP/jeter.vcf.gz

    mv TMP/jeter.vcf.gz ./${prefix}.vcf.gz
    mv TMP/jeter.vcf.gz.tbi ./${prefix}.vcf.gz.tbi
    """
}

workflow STRUCTURAL_VARIANTS {
    take:
        meta
        fasta
        fai
        dict
        gff3
        pedigree
        vcfsv
    main:

        TRUVARI(
	        [id:"truvari"],
            fasta,
            fai,
            dict,
            vcfsv
        )
        ANNOTATE_SV(
            meta,
            fasta,
            fai,
            dict,
            gff3,
            TRUVARI.out.vcf
        )
        
        NO_CALL_TO_HOM_REF(ANNOTATE_SV.out.vcf)
        BCTOOLS_MENDELIAN2(fai, pedigree, NO_CALL_TO_HOM_REF.out.vcf)

}

workflow HET_COMPOSITE {
    take:
        meta
        fasta
        fai
        dict
        pedigree
        vcfs
    main:
        FILTER_FOR_HET_COMPOSITE(vcfs)
}

process FILTER_FOR_HET_COMPOSITE {
    tag "${meta.id}"
    label "process_single"
    conda "${moduleDir}/../../conda/bioinfo.01.yml"
    input:
        tuple val(meta),path(vcf),path(idx)
    script:
        def args = meta.contig?"--regions \"${meta.contig}\"":""
        def prefix = task.ext.prefix?:meta.id+(meta.contig?"."+meta.contig:"")
        """

| gunzip -c | java -jar ~/jvarkit.jar vcffilterjdk -e 'return getVepPredictionParser().getPredictions(variant).stream().filter(P->P.get("Allele").matches("[ATGC]+)).map(P->P.get("gnomADg_AF_nfe")).filter(it->it!=null && !it.trim().isEmpty()).mapToDouble(it->Double.parseDouble(it)).filter(it->it<=1E-100).count()==1L;' | java -jar ~/jvarkit.jar vcf2table --hide "GENOTYPES"



        """
    }

