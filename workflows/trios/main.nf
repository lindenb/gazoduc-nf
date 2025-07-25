/*

Copyright (c) 2025 Pierre Lindenbaum

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

include {VEP                                      } from '../../subworkflows/annotation/vep/main.nf'
include {DOWNLOAD_GTF_OR_GFF3 as DOWNLOAD_GFF3    } from '../../modules/gtf/download/main.nf'
include {DOWNLOAD_GTF_OR_GFF3 as DOWNLOAD_GTF     } from '../../modules/gtf/download/main.nf'
include {BCFTOOLS_BCSQ                            } from '../../modules/bcftools/bcsq/main.nf'
include {BCTOOLS_MENDELIAN2                       } from '../../modules/bcftools/mendelian2/main.nf'
include {HET_COMPOSITE  as APPPLY_HET_COMPOSITE   } from '../../modules/jvarkit/hetcomposite/main.nf'
include {BCFTOOL_CONCAT as CONCAT1                } from '../../modules/bcftools/concat/main.nf'
include {WORKFLOW_SNV                             } from './sub.snv.nf'
include {WORKFLOW_SV                              } from './sub.sv.nf'
include {SOMALIER_BAMS                            } from '../../subworkflows/somalier/bams/main.nf'
include {MULTIQC                                  } from '../../modules/multiqc'
include {COMPILE_VERSIONS                         } from '../../modules/versions/main.nf'

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
        def ref_hash = [
            id: file(params.fasta).simpleName,
            ucsc_name: (params.ucsc_name?:"undefined")
            ]
        def fasta  = [ref_hash, file(params.fasta) ]
        def fai    = [ref_hash, file(params.fai) ]
        def dict   = [ref_hash, file(params.dict) ]
        def bed    = params.bed!=null ? [ref_hash,file(params.bed)] : [ref_hash, [] ]
        def pedigree = [[id:"pedigree"],file(params.pedigree)]
        
        def versions = Channel.empty()
        to_multiqc = Channel.empty()

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

        bams = Channel.empty()
        triosbams_ch = Channel.empty()
        
        if(params.bams) {
            bams = Channel.fromPath(params.bams)
                .splitCsv(header:true,sep:',')
                .map{assertKeyMatchRegex(it,"sample","^[A-Za-z_0-9\\.\\-]+\$")}
                .map{assertKeyMatchRegex(it,"bam","^\\S+\\.(bam|cram)\$")}
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
                .map{[[id:it.sample],file(it.bam),file(it.bai),file(it.fasta),file(it.fai),file(it.dict)]}
            

            SOMALIER_BAMS(
                [id:"somalier"],
                fasta,
                fai,
                dict,
		        triosbams_ch, // sample,bam,bai
		        pedigree, // pedigree for somalier
		        [[id:"no_sites"],[]]
                )
            versions = versions.mix(SOMALIER_BAMS.out.versions)

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

        DOWNLOAD_GFF3(fasta,fai,dict)
        versions = versions.mix(DOWNLOAD_GFF3.out.versions)

        DOWNLOAD_GTF(fasta,fai,dict)
        versions = versions.mix(DOWNLOAD_GTF.out.versions)
    /*
        STRUCTURAL_VARIANTS(
            [id:"sv"],
            fasta,
            fai,
            dict,
            DOWNLOAD_GTF.out.output,
            pedigree,
            vcfs.sv.take(10)//TODO
        )*/



        WORKFLOW_SNV(
            [id:"snv"],
            fasta,
            fai,
            dict,
            bed,
            pedigree,
            [ref_hash,file(params.gnomad),file(params.gnomad+".tbi")],
            DOWNLOAD_GFF3.out.output,
            DOWNLOAD_GTF.out.output,
            triosbams_ch,
            vcfs.snv
            )
        versions = versions.mix(WORKFLOW_SNV.out.versions)

        WORKFLOW_SV(
            [id:"snv"],
            fasta,
            fai,
            dict,
            bed,
            pedigree,
            DOWNLOAD_GFF3.out.output,
            DOWNLOAD_GTF.out.output,
            triosbams_ch,
            vcfs.sv
            )
            
       versions = versions.mix(WORKFLOW_SV.out.versions)


	    COMPILE_VERSIONS(versions.collect())
	    to_multiqc = to_multiqc.mix(COMPILE_VERSIONS.out.multiqc)

       MULTIQC(to_multiqc.collect().map{[[id:"trios"],it]})
}


process FILTER_MENDELIAN {
    tag "${meta.id?:vcf.name}"
    label "process_single"
    conda "${moduleDir}/../../conda/bioinfo.01.yml"
    afterScript "rm -rf TMP"
    input:
        tuple val(meta1 ),path(pedigree)
        tuple val(meta ),path(vcf),path(idx)
    output:
        tuple val(meta),path("*.bcf"),path("*.csi"),emit:vcf
    script:
        def prefix = task.ext.prefix?:vcf.baseName+".merr"
    """
    mkdir -p TMP
cat << EOF > TMP/jeter.code
final Set<String> controls = new HashSet<>(Arrays.asList(
EOF

awk '(\$3=="0" && \$4=="0") {print \$2;}' "${pedigree}" | sort| uniq > TMP/jeter.controls

comm -13 <(cut -f2  "${pedigree}" | sort | uniq) \\
     <(bcftools query -l "${vcf}"| sort | uniq)  >> TMP/jeter.controls


 awk '{printf("\"%s\"\\n",\$2);}' TMP/jeter.controls | paste -s -d, >> TMP/jeter.code

cat << EOF >> TMP/jeter.code
));
return variant
    .getGenotypes()
    .stream()
    .filter(G->controls.contains(G.getSampleName()))
    .noneMatch(G->G.hasAltAllele())
    ;
EOF


    bcftools view -i 'INFO/loConfDeNovo!="" || INFO/hiConfDeNovo!="" || INFO/MERR>0' '${vcf}' |\\
        jvarkit -XX:-UsePerfData -Djava.io.tmpdir=TMP vcffilterjdk --filter 'CONTROL_WITH_ALT' -f TMP/jeter.code |\\
        bcftools view -O b -o ${prefix}.bcf
    bcftools index --threads "${task.cpus}"  -f "${prefix}.bcf"
    """
}

process MERGE_AND_FILTER {
    tag "${meta.id} ${meta.contig?:""}"
    label "process_single"
    conda "${moduleDir}/../../conda/bioinfo.01.yml"
    afterScript "rm -rf TMP"
    input:
        tuple val(met1),path(pedigree)
        tuple val(meta),path("VCFS/*")
    output:
        tuple val(meta),path("*.vcf.gz"),path("*.tbi"),emit:vcf
        tuple val(meta),path("*.table.txt")
        tuple val(meta),path("*.genes.tsv")
        tuple val(meta),path("*.bed"),emit:bed //used for igv reports
    script:
        def prefix = task.ext.prefix?:"snv"
    """
    mkdir -p TMP
    find VCFS/  -name "*.vcf.gz" -o -name "*.bcf" | sort > TMP/jeter.list
    
    bcftools concat --allow-overlaps --file-list TMP/jeter.list -O z -o TMP/jeter.vcf.gz
    bcftools index -t -f TMP/jeter.vcf.gz

    bcftools view --apply-filters '.,PASS' TMP/jeter.vcf.gz |\\
        jvarkit -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP vcf2table \\
            --hide 'HOM_REF,NO_CALL' --pedigree ${pedigree} > TMP/jeter.table.txt

    bcftools view --apply-filters '.,PASS' TMP/jeter.vcf.gz |\\
        jvarkit -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP groupbygene > TMP/jeter.genes.tsv

    bcftools view --apply-filters '.,PASS' TMP/jeter.vcf.gz |\\
        jvarkit -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP bioalcidaejdk \\
            -e 'stream().forEach{
                final Set<String> set = new HashSet<>();
                for(int side=0;side<2;side++) {
                    final String tag=(i==0?"hiConfDeNovo":"loConfDeNovo");
                    if(!variant.hasAttribute(tag)) continue; 
                    set.addAll(variant.getAttributeAsStringList(tag,""));
                    }
                for(String s : set) {
                    println(variant.getContig()+"\t"+(variant.getStart()-1)+"\t"+variant.getEnd()+"\t"+s);
                    }
                }' |\\
        LC_ALL=C sort -T TMP -k1,1 -k2,2n |\\
        uniq > TMP/jeter.bed

    mv TMP/jeter.vcf.gz ./${prefix}.vcf.gz
    mv TMP/jeter.vcf.gz.tbi ./${prefix}.vcf.gz.tbi
    mv TMP/jeter.table.txt ./${prefix}.table.txt
    mv TMP/jeter.genes.tsv ./${prefix}.genes.tsv
    mv TMP/jeter.bed  ./${prefix}.bed
    """
}

/********************************************************************************************
 * STRUCTURAL_VARIANTS
 */
workflow STRUCTURAL_VARIANTS {
    take:
        meta
        fasta
        fai
        dict
        gtf
        pedigree
        vcfsv
    main:
        versions = Channel.empty()
      
      
      
        //BCTOOLS_MENDELIAN2(fai, pedigree, NO_CALL_TO_HOM_REF.out.vcf)
	vcf = Channel.empty()
    emit:
	vcf
        //vcf = BCTOOLS_MENDELIAN2.out.vcf
        versions
}

/********************************************************************************************
 * HET_COMPOSITE
 */
workflow HET_COMPOSITE {
    take:
        meta
        fasta
        fai
        dict
        pedigree
        vcfs
    main:
        versions = Channel.empty()
        /** merge VCFS */
        CONCAT1(
            vcfs
                .map{[it[1],it[2]]}
                .collect()
                .map{[meta,it.flatten()]},
            [[:],[]]//NO BED
            )
        versions = versions.mix(CONCAT1.out.versions)

        FILTER_FOR_HET_COMPOSITE(CONCAT1.out.vcf)
        versions = versions.mix(FILTER_FOR_HET_COMPOSITE.out.versions)
        
        APPPLY_HET_COMPOSITE(
                fasta,
                fai,
                dict,
                pedigree,
                FILTER_FOR_HET_COMPOSITE.out.vcf
                )
        versions = versions.mix(APPPLY_HET_COMPOSITE.out.versions)
        
    emit:
        vcf = APPPLY_HET_COMPOSITE.out.vcf
        versions
}

process FILTER_FOR_HET_COMPOSITE {
tag "${meta.id}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
    tuple val(meta),path(vcf),path(vcfix)
output:
    tuple val(meta),path("*.bcf"),path("*.csi"),emit:vcf
    path("versions.yml"),emit:versions

script:
    def prefix = task.ext.prefix?:meta.id+".filterhetcomposite"
"""
mkdir -p TMP
bcftools view ${vcf} |\\
    bcftools view --apply-filters '.,PASS' --write-index -O b -o ${prefix}.bcf

cat << END_VERSIONS > versions.yml
"${task.process}":
    bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
}

