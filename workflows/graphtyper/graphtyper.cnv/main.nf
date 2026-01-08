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
nextflow.enable.dsl=2

include {GRAPHTYPER_GENOTYPE_SV           } from '../../../modules/graphtyper/genotype_sv'
include {runOnComplete;dumpParams         } from '../../../modules/utils/functions.nf'
include {removeCommonSuffixes             } from '../../../modules/utils/functions.nf'
include {MULTIQC                          } from '../../../modules/multiqc'
include {COMPILE_VERSIONS                 } from '../../../modules/versions/main.nf'
include {JVARKIT_VCF_SET_DICTIONARY       } from '../../../modules/jvarkit/vcfsetdict'
include {JVARKIT_VCF_CLUSTER              } from '../../../modules/jvarkit/vcfcluster'
include {BCFTOOLS_CONCAT                  } from '../../../modules/bcftools/concat'
include {JVARKIT_BAM_RENAME_CONTIGS       } from '../../../modules/jvarkit/bamrenamechr'
include {PLOT_COVERAGE_01                 } from '../../../subworkflows/plotdepth'
include { PREPARE_ONE_REFERENCE           } from '../../../subworkflows/samtools/prepare.one.ref'
include {READ_SAMPLESHEET as READ_BAMS    } from '../../../subworkflows/nf/read_samplesheet'
include { META_TO_BAMS                    } from '../../../subworkflows/samtools/meta2bams1'
include { PREPARE_USER_BED                } from '../../../subworkflows/bedtools/prepare.user.bed'

if(params.help) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}



workflow {
	versions = Channel.empty()
	multiqc = Channel.empty()


    if(params.fasta==null) {
        throw new IllegalArgumentException("undefined --fasta");
        }
	if(params.samplesheet==null) {
		throw new IllegalArgumentException("undefined --samplesheet");
		}
    if(params.vcf==null) {
		throw new IllegalArgumentException("undefined --vcf");
		}
	multiqc = Channel.empty()
	versions = Channel.empty()
	def metadata = [id:"graphtyper.cnv"]
		

	PREPARE_ONE_REFERENCE(
		metadata,
		Channel.fromPath(params.fasta).map{f->[[id:f.baseName],f]}
		)
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)
	


	//def gtf  =    [ metadata, file(params.gtf),file(params.gtf +".tbi") ]
	//def pedigree = [ metadata, []]

	READ_BAMS(
		[arg_name:"samplesheet"],
		params.samplesheet
		)
	versions = versions.mix(READ_BAMS.out.versions)
	

	META_TO_BAMS(
		metadata ,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		READ_BAMS.out.samplesheet
		)
	versions = versions.mix(META_TO_BAMS.out.versions)



    

	def vcf = Channel.of(params.vcf)
        .map{vcf-> [
            [id:removeCommonSuffixes(file(vcf).name)],
            file(vcf)
            ]}


	JVARKIT_VCF_CLUSTER( vcf )
	versions = versions.mix(JVARKIT_VCF_CLUSTER.out.versions)

    vcfs = JVARKIT_VCF_CLUSTER.out.vcf
            .map{meta,f->(f instanceof List?f:[f])}
            .flatMap()
            .map{f->[f.name+".tbi",f]}
    tbis = JVARKIT_VCF_CLUSTER.out.tbi
            .map{meta,f->(f instanceof List?f:[f])}
            .flatMap()
            .map{f->[f.name,f]}
    vcfs = vcfs.join(tbis)
        .map{_meta,f,tbi->[
            [id:removeCommonSuffixes(f.name)],
            f,
            tbi
            ]}




    grouped_bams = META_TO_BAMS.out.bams.map{meta,bam,bai->[bam,bai]}
        .flatMap()
        .collect()
        .map{[[id:"bams"],it.sort()]}


    GRAPHTYPER_GENOTYPE_SV(
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		grouped_bams,
		vcfs
		)

    /*
    vcfs = JVARKIT_VCF_CLUSTER.out.vcf
        .map{meta,files->[meta,(files instanceof List?files:[files])]}
        .flatMap()
        .map{vcf->[[id:removeCommonSuffixes(vcf.name)],vcf]}

	 = bams2_ch.
   
	
	versions = versions.mix(GRAPHTYPER_GENOTYPE_SV.out.versions)
    */

    
	/*

	BCFTOOLS_CONCAT(
       	GRAPHTYPER_GENOTYPE_SV.out.vcf
            .map{[it[1],it[2]]}
            .collect()
            .map{[[id:"gtyper"],it]},
        [[id:"nobed"],[]]
        )
     versions =  versions.mix(BCFTOOLS_CONCAT.out.versions)


	JVARKIT_VCF_SET_DICTIONARY(
        dict,
        BCFTOOLS_CONCAT.out.vcf
        )
    versions = versions.mix(JVARKIT_VCF_SET_DICTIONARY.out.versions)


    if(params.pedigree!=null) {
        GTYPER_DENOVO(
            fasta,
            fai,
            dict,
            pedigree,
            JVARKIT_VCF_SET_DICTIONARY.out.vcf
            )
        versions = versions.mix(GTYPER_DENOVO.out.versions)

        EXTRACT_DENOVO( GTYPER_DENOVO.out.vcf);
        versions = versions.mix(EXTRACT_DENOVO.out.versions)


	chcov = EXTRACT_DENOVO.out.bed.map{it[1]}.splitCsv(sep:'\t',header:true)
		.map{[
			contig:it.chrom,
			start: it.start,
			end: it.end,
			title: it.title,
			id : it.chrom+"_"+it.start+"_"+it.end+"."+ it.toString().md5().substring(0,7)
			]}
		.combine(Channel.of(fasta).map{it[1]})
		.combine(Channel.of(fai).map{it[1]})
		.combine(Channel.of(dict).map{it[1]})
		.combine(grouped_bams.map{[it[1]]}) */ /* convert to array */
    /*
	PLOT_COVERAGE_01(
        [id:"graphtyper"],
        fasta,
        fai,
        dict,
        gtf,
        EXTRACT_DENOVO.out.bed,
        bams2_ch
        )
    versions = versions.mix(PLOT_COVERAGE_01.out.versions)
    */
    }
	

    

    /*
	COMPILE_VERSIONS(versions.collect())
       multiqc = multiqc.mix(COMPILE_VERSIONS.out.multiqc)

    MULTIQC(multiqc.collect().map{[[id:"gtyper_sv"],it]})
	*/
	

runOnComplete(workflow);

process GTYPER_DENOVO {
    label "process_single"
    tag "${meta.id?:""} ${vcf.name}"
    afterScript "rm -rf TMP"
    conda "${moduleDir}/../../../conda/bioinfo.01.yml"
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
	"""
	hostname 1>&2
	mkdir -p TMP

cat << __EOF__ > TMP/jeter.code

private boolean isFiltered(final Genotype g) {
	if(!g.isFiltered()) return false;
	if(g.getFilters().equals("PASS")) return false;
	return true;
 	}

private boolean acceptControl(final VariantContext vc,String sm) {
        final Genotype g = vc.getGenotype(sm);
        if(g!=null && g.hasAltAllele() ) return false;
        return true;
        }

private boolean acceptTrio(final VariantContext vc,String cm,String fm,String mm) {
    final Genotype c = vc.getGenotype(cm);
    if(c==null) return false;
    if(isFiltered(c)) return false;
    final Genotype m = vc.getGenotype(mm);
    final Genotype f = vc.getGenotype(fm);
    // child must be HET
    if(!c.hasAltAllele()) return false;
    // keep de novo
    if(f.hasAltAllele() || m.hasAltAllele()) return false;
    if(isFiltered(f)) return false;
    if(isFiltered(m)) return false;
    return true;
    }

public Object apply(final VariantContext variant) {
if(variant.isFiltered()) return false;
final String svType = variant.getAttributeAsString("SVTYPE","");
if(svType.equals("BND")) return false;
__EOF__



    # convert pedigree if no 6th column
    awk '{S=\$6 ; if(NF==5 || S=="") { if(\$3!="0" && \$4!="0") {S="case";} else {S="control"} }  printf("%s\t%s\t%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,\$4,\$5,S);}' "${pedigree}" > TMP/pedigree.tsv
    

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


    # GATK possibledenovo
    awk -f "${moduleDir}/../../../modules/gatk/possibledenovo/pedigree4gatk.awk" "${pedigree}" > TMP/jeter.ped

    gatk --java-options "-XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" VariantAnnotator \\
        -R "${fasta}" \\
        --annotation PossibleDeNovo \\
        --pedigree  TMP/jeter.ped \\
        -V TMP/jeter.vcf.gz \\
        -O TMP/jeter2.vcf.gz




    mv  TMP/jeter2.vcf.gz "${prefix}.vcf.gz"
    mv  TMP/jeter2.vcf.gz.tbi "${prefix}.vcf.gz.tbi"


cat << EOF > versions.yml
"${task.process}":
	jvarkit: todo
	gatk: todo
EOF
"""
}






process EXTRACT_DENOVO {
    label "process_single"
    tag "${meta.id?:""}"
    afterScript "rm -rf TMP"
    conda "${moduleDir}/../../../conda/bioinfo.01.yml"
    input:
	    tuple val(meta),path(vcf),path(vcfidx)
    output:
        tuple val(meta),path("*.bed"),emit:bed
        path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:"\${MD5}.denovo"
"""
mkdir -p TMP

echo -e 'chrom\tstart\tend\ttitle' >> TMP/jeter.bed

bcftools query -f '%CHROM\t%POS0\t%END\t%SVTYPE %hiConfDeNovo %loConfDeNovo\\n' "${vcf}" |\\
    tr "," " " | tr -s " " |\\
   sort -t '\t' -k1,1 -k2,2n  -k3,3n -T TMP --unique   >> TMP/jeter.bed

MD5=`cat TMP/jeter.bed | md5sum | cut -d ' ' -f1`

mv TMP/jeter.bed "${prefix}.bed"

cat << EOF > versions.yml
"${task.process}":
        jvarkit: todo
        gatk: todo
EOF
"""
}
