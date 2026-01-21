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
process COMPILE_VERSIONS {
label "process_single"
afterScript "rm -rf TMP"
input:
    path("versions??.yml")
output:
    path("*_mqc.tsv"),emit:multiqc
script:
    def prefix=task.ext.prefix?:"versions"
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
    awk 'BEGIN{printf("process\tversions\\n");} {print}' > "${prefix}_mqc.tsv"
"""

stub:
"""
touch versions_mqc.tsv
"""
}
