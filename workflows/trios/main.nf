import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

include {CLINVAR                } from '../../subworkflows/annotation/clinvar/main.nf'
include {ALPHAMISSENSE          } from '../../subworkflows/annotation/alphamissense/main.nf'
include {VEP                    } from '../../subworkflows/annotation/vep/main.nf'
include {TRIOS                  } from '../../subworkflows/trios/main.nf'
include {TRUVARI                } from '../../subworkflows/truvari/main.nf'
include {SNPEFF                 } from '../../subworkflows/snpeff/main.nf'

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

//        MAKE_PED([[id:"pedigree"],file(params.pedigree)])

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

        TRUVARI(
	        [id:"truvari"],
            fasta,
            fai,
            dict,
            vcfs.sv
        )


        vcf2bed = VCF_TO_CONTIGS(vcfs.snv) 
        vcf_snv = vcf2bed.output
            .splitText(elem:1)
            .map{[it[0],it[1].trim(),it[2],it[3]]}
            .map{[it[0].plus(contig:it[1]),it[2],it[3]]}.view()

        
        SUBSET_VCF(vcf_snv.view{"SUBSET $it \n"})
        vcf_snv = SUBSET_VCF.out.vcf

        VEP([id:"vep"],fasta,fai,dict,vcf_snv)
        vcf_snv = SUBSET_VCF.out.vcf

        CLINVAR([id:"clinvar"],fasta,fai,dict,NO_BED,VEP.out.vcf)
        vcfs = CLINVAR.out.vcf
 
        ALPHAMISSENSE([id:"alphamissense"],fasta,fai,dict,NO_BED,CLINVAR.out.vcf)
        vcfs = ALPHAMISSENSE.out.vcf

        SNPEFF([id:"snpeff"],fasta,fai,dict,ALPHAMISSENSE.out.vcf)
        vcfs = SNPEFF.out.vcf
        
/*
        TRIO([id:"trio"],fasta,fai,dict,SNPEFF.out.vcf)
        vcfs = VEP.out.vcf
	*/
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
        tuple val(meta),path(vcf),path(idx)
    output:
        tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
    script:
        def args = meta.contig?"--regions \"${meta.contig}\"":""
        def prefix = task.ext.prefix?:meta.id+(meta.contig?"."+meta.contig:"")
        """
        bcftools view --threads "${task.cpus}" -O z ${args} -o ${prefix}.vcf.gz "${vcf}" 
        bcftools index --threads "${task.cpus}" -t -f ${prefix}.vcf.gz
        """
}

process MAKE_PED {
 tag "${meta.id}"
    label "process_single"
    input:
        tuple val(meta),path(ped)
    output:
        tuple val(meta),path("pedigree.tsv"),emit:output
    script:
    """
    """
}

wokflow HET_COMPOSITE {
    take:
        meta
        fasta
        fai
        dict
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
        jvarkit vcffilterjdk -e 'getTools(). ().'
        """
    }

