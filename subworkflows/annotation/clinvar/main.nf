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

include {k1_signature} from '../../../modules/utils/k1.nf'


workflow CLINVAR {
	take:
                meta
		fasta
		fai
		dict
                bed
		vcfs /* meta, vcf,vcf_index */
	main:
                versions = Channel.empty()
                DOWNLOAD(
                        fasta,
                        fai,
                        dict,
                        bed
                        )
                versions = versions.mix(DOWNLOAD.out.versions)

		ANNOTATE(DOWNLOAD.out.vcf, vcfs)
                versions = versions.mix(ANNOTATE.out.versions)
	emit:
		vcf = ANNOTATE.out.vcf
                versions
	}


process DOWNLOAD {
afterScript "rm -rf TMP"
label "process_quick"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta4),path(bed)
output:
    tuple val(meta1),path("*.bcf"),path("*.bcf.csi"),emit:vcf
    path("versions.yml"),emit:versions
script:
    def k1= k1_signature();
    def prefix = "clinvar"
    def base = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar"
    def args = bed?"--regions-overlap 1 --targets-file TMP/nochr.bed":""

"""
set -o pipefail
set -x
mkdir -p TMP


if ${bed?true:false}
then
       sed 's/^chr//' '${bed}' >  TMP/nochr.bed
fi

cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg38}\t${base}/vcf_GRCh38/clinvar.vcf.gz
1:${k1.hg19}\t${base}/vcf_GRCh37/clinvar.vcf.gz
EOF


awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sed 's/^chr//' | sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv
join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv | sort | uniq > TMP/jeter.url

test -s TMP/jeter.url

wget -O - `cat TMP/jeter.url` |\\
        bcftools view -O v ${args}  |\\
        jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  vcfsetdict -R "${fasta}"  -n SKIP |\\
        bcftools view -O b -o TMP/jeter.bcf

bcftools view --header-only TMP/jeter.bcf | grep "^##INFO" | cut -d '=' -f 3 | cut -d, -f 1| grep -v "#"  |\\
        awk '{printf("INFO/%s\tCLINVAR_%s\\n",\$1,\$1);}' > TMP/rename.tsv

bcftools annotate \\
    --threads ${task.cpus} \\
    --rename-annots TMP/rename.tsv \\
    -O b \\
    -o TMP/${prefix}.bcf \\
    TMP/jeter.bcf

bcftools index  \\
    --threads ${task.cpus} \\
     TMP/${prefix}.bcf

mv TMP/${prefix}.bcf ./
mv TMP/${prefix}.bcf.csi ./


cat << EOF > versions.yml
${task.process}:
    bcftools: TODO
EOF
"""
}

process ANNOTATE {
tag "${meta.id}"
label "process_quick"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
        tuple val(meta1),path(clinvar),path(clinvaridx)
        tuple val(meta ),path(vcf),path(tbi)
output:
        tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
        path("versions.yml"),emit:versions
script:
        def prefix = task.ext.prefix?:vcf.baseName+".clinvar"
"""
mkdir -p TMP

bcftools annotate \\
        --threads ${task.cpus} \\
        -a "${clinvar}" \\
        -c "CLINVAR_AF_ESP,CLINVAR_AF_EXAC,CLINVAR_AF_TGP,CLINVAR_ALLELEID,CLINVAR_CLNDN,CLINVAR_CLNDNINCL,CLINVAR_CLNDISDB,CLINVAR_CLNDISDBINCL,CLINVAR_CLNHGVS,CLINVAR_CLNREVSTAT,CLINVAR_CLNSIG,CLINVAR_CLNSIGCONF,CLINVAR_CLNSIGINCL,CLINVAR_CLNVC,CLINVAR_CLNVCSO,CLINVAR_CLNVI,CLINVAR_DBVARID,CLINVAR_GENEINFO,CLINVAR_MC,CLINVAR_ONCDN,CLINVAR_ONCDNINCL,CLINVAR_ONCDISDB,CLINVAR_ONCDISDBINCL,CLINVAR_ONC,CLINVAR_ONCINCL,CLINVAR_ONCREVSTAT,CLINVAR_ONCCONF,CLINVAR_ORIGIN,CLINVAR_RS,CLINVAR_SCIDN,CLINVAR_SCIDNINCL,CLINVAR_SCIDISDB,CLINVAR_SCIDISDBINCL,CLINVAR_SCIREVSTAT,CLINVAR_SCI,CLINVAR_SCIINCL" \\
        -O z \\
        -o TMP/jeter.vcf.gz \\
        ${vcf}

bcftools index --threads ${task.cpus} -t -f "TMP/jeter.vcf.gz"

mv TMP/jeter.vcf.gz "${prefix}.vcf.gz"
mv TMP/jeter.vcf.gz.tbi "${prefix}.vcf.gz.tbi"

cat << EOF > versions.yml
${task.process}:
	bcftools: TODO
EOF
"""
}
