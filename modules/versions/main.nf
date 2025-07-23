process COMPILE_VERSIONS {
label "process_single"
afterScript "rm -rf TMP"
input:
    path("versions??.yml")
output:
    path("versions_mqc.tsv"),emit:multiqc
script:
"""
mkdir -p TMP
cat << 'EOF' > jeter.awk
/^#/ {
    next;
    }

/^[\t ]*\$/ {
    next;
    }

/^[\t ]+/ {
    L=\$0;
    gsub(/^[\t ]*/,"",L);
    gsub(/[\t ]*\$/,"",L);
    gsub(/[\t ]/," ",L);
    printf("%s\t%s\t%s\\n",TIME,STEP,L);
    next;
    }

/^[^ \t]/ {
    STEP=\$0;
    gsub(/^[:\\"]+/,"",STEP);
    gsub(/[\t :\\"]+\$/,"",STEP);
    }
 
EOF

find . -name "versions*.yml" -printf  "%T@\t%p\\n" |\\
while read TIME YAML
do
    awk -v TIME=\${TIME} -f jeter.awk "\${YAML}" 
done |\\
    sort -t '\t' -T TMP -k2,2 -k3,3 --unique |\\
    sort -t '\t' -T TMP -k1.1g |\\
    cut -f2,3 |\\
    awk 'BEGIN{printf("process\tversions\\n");} {print}' > versions_mqc.tsv
"""
}
