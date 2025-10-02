
process REDUCE_BED {
tag "${bed.name}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta),path(bed)
output:
    tuple val(meta),path("BEDS/*"),optional:true,emit:bed
    path("versions.yml"),emit:versions
script:
    def njobs  = task.ext.njobs?:10
    def prefix = task.ext.prefix?:"${meta.id}"
"""
mkdir -p BEDS

# not empty please
test -s "${bed}"

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
        split -a 9 --additional-suffix=.bed --lines=1 - BEDS/${prefix}.${bed.name.md5()}

fi

touch versions.yml
"""
}
