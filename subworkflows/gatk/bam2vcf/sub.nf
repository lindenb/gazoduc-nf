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
include {GATK_BAM2VCF   } from '../../../modules/gatk/bam2vcf'


String toFilePath(Object bed) {
	def dir2 = bed.toRealPath().toString().md5();
	def f= ""+workflow.workDir + "/FLAG/" + dir2.substring(0,2) + "/" + dir2.substring(2,4);
    return f;
	}


String toFailurePath(Object bed) {
	return toFilePath(bed)+".failed";
	}
String toSuccessPath(Object bed) {
	return toFilePath(bed)+".ok";
	}


workflow DIVIDE_AND_CONQUER {
take:
    meta
    level
    fasta
    fai
    dict
    dbsnp
    references
    bed // [meta,bed]
    bams // [meta, [bams and bai] ]
main:
    versions = Channel.empty()
    ok_ch = Channel.empty()
    bed = bed.map{[it[0].plus("level":level),it[1]]}

    bed.branch {B->
        known_for_success : file(toSuccessPath(B[1])).exists()
        known_to_fail: file(toFailurePath(B[1])).exists()
        todo: true
        }.set{branch1}

   //branch1.known_to_fail.view{"Known to fail: ${it}"}
   //branch1.todo.view{"TODO to fail: ${it}"}

    GATK_BAM2VCF(
        fasta,
        fai,
        dict,
        dbsnp,
        references,
        bams.combine(branch1.todo)
        .map{[
            it[2]/* bed.meta */,
            it[1]/*bam and bai */,
            it[3]/*bed */
            ]}
        )
    versions = versions.mix(GATK_BAM2VCF.out.versions)


    branch1.todo /* meta,bed */
        .mix(GATK_BAM2VCF.out.vcf.map{[it[0],it[3]]}) /* meta bed */
        .groupTuple()
        .branch{v->
            success: v[1].size()==2
            failure: v[1].size()==1
            other: true
            }.set{branch2}
    
    branch2.other.map{throw new IllegalArgumentException("${it}");}



    TOUCH_FAILURE(branch2.failure)

    SPLITBED(
        fasta,
        fai,
        dict,
        level,
        TOUCH_FAILURE.out.bed.mix(branch1.known_to_fail)
        )
    
    todo_bed = SPLITBED.out.bed
        .map{it[1]}
        .map{it instanceof List?it:[it]}
        .flatMap()
        .map{[[id:it.toRealPath().toString().md5()],it]}      

    TOUCH_SUCCESS(GATK_BAM2VCF.out.vcf)

    success_vcf = branch1.known_for_success
        .map{file(toSuccessPath(it[1]))}
        .splitCsv(header:false,sep:',')
        .map{[[id:it[0]],file(it[1]),file(it[2])]}

    sucecss_vcf = success_vcf.mix(TOUCH_SUCCESS.out.vcf)

emit:
    versions
    bed = todo_bed
    vcf = success_vcf
}

process TOUCH_FAILURE {
tag "${bed.name}"
executor "local"
input:
    tuple val(meta),path(bed)
output:
    tuple val(meta),path(bed),emit:bed
script:
    def f = file(toFailurePath(bed));
"""
mkdir -p "${f.parent}"
echo "${bed.toRealPath()}" > ${f}
"""
}

process TOUCH_SUCCESS {
tag "${bed.name}"
executor "local"
input:
    tuple val(meta),path(vcf),path(tbi),path(bed)
output:
    tuple val(meta),path(vcf),path(tbi),emit:vcf
script:
    def f = file(toSuccessPath(bed));
"""
mkdir -p "${f.parent}"
echo '${meta.id},${vcf.toRealPath()},${tbi.toRealPath()}' > ${f}
"""
}



process SPLITBED {
tag "${bed.name} Level ${level}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    val(level)
    tuple val(meta),path(bed)
output:
    tuple val(meta),path("BEDS/*"),optional:true,emit:bed
    path("versions.yml"),emit:versions
script:
    def njobs  = task.ext.njobs?:"10"
"""
mkdir -p BEDS
if [[ \$(wc -l < "${bed}") -gt 1  ]]
then
    jvarkit  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData bedcluster \\
        -R ${fasta} \\
        --jobs ${njobs} \\
        -o BEDS
else
    awk -F '\t' '{
        B=int(\$2);
        E=int(\$3);
        N=E-B;
        if(N<=1) next;
        DX=int(N/${njobs});
        if(DX<1) DX=1;
        while(B < E) {
            E2 = B+DX;
            if(E2 > E) E2 = E;
            printf("%s\t%d\t%d\\n",\$1,B,E2);
            B+=DX;
            }
        }' ${bed} |\\
        split -a 9 --additional-suffix=.bed --lines=1 - BEDS/LEVEL.${level}.${bed.name.md5()}

fi

touch versions.yml
"""
}