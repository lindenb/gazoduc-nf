
process SPLIT_BED {
tag "${bed.name}"
executor "local"
maxForks 10
input:
        tuple path(bed),path(vcfs)
output:
        tuple path("BEDS/*.bed"),path(vcfs),emit:output
script:
        def f = bed.getBaseName();
        def size=task.ext.size?:100000;
"""
mkdir -p BEDS
awk -F '\t' 'BEGIN{N=1;T=0.0;f=sprintf("BEDS/${f}.%d.N${size}.bed",N);} {print \$0 >> f; T+=int(\$3)-int(\$2); if(T>=${size}) {close(f);T=0.0;N++;f=sprintf("BEDS/${f}.%d.N${size}.bed",N);}}' '${bed}'
"""
}
