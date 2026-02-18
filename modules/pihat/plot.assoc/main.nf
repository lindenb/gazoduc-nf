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
process PLOT_ASSOC {
tag "${assoc.name}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
	tuple val(meta ),path(assoc)
output:
    tuple val(meta ),path("*.png"),emit:png
    path("versions.yml"),emit:versions
script:
"""
hostname 1>&2
mkdir -p TMP
set -x

JD1=`which jvarkit`
echo "JD1=\${JD1}" 1>&2
# directory of jvarkit
JD2=`dirname "\${JD1}"`
# find the jar itself
JVARKIT_JAR=`find "\${JD2}/../.." -type f -name "jvarkit.jar" -print -quit`
JVARKIT_DIST=`dirname "\${JVARKIT_JAR}"`


cat "${moduleDir}/Minikit.java" |\\
    sed -e 's/__INPUT__/${assoc}/'  -e 's/__MIN_P_VALUE__/1e-8/' -e 's/__FASTA__/${fasta}/'  > TMP/Minikit.java

javac -d TMP -cp \${JVARKIT_DIST}/jvarkit.jar TMP/Minikit.java
java -Djava.awt.headless=true -cp \${JVARKIT_DIST}/jvarkit.jar:TMP Minikit

#
# Floriane suggested to also add the number of variants per chromosome
#
awk '(\$2!="SNP") {print \$2}' '${assoc}' |\\
	cut -d ':' -f 1 | sort -T TMP | uniq -c |\\
	awk '{printf("%s\t%s\\n",\$2,\$1);}' |\\
	sort -T TMP -t '\t' -k1,1 |\\
	join -t '\t' -1 1 -2 1 - <(sort -t '\t' -k1,1 '${fai}' ) |\\
	awk 'BEGIN{printf("<table class=\\"table\\"><thead><caption>Number of variants per contig</caption><tr><th>CHROM</th><th>Length</th><th>VARIANTS</th><th>Variants per base</th></tr></thead><tbody>\\n");} {printf("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>\\n",\$1,\$3,\$2,\$2/\$3);} END {printf("</tbody></table>\\n");}' >> TMP/jeter.html

mv TMP/jeter.png "${assoc.baseName}.png"

cat << EOF > versions.yml
${task.process}:
    R: todo
EOF
"""
}
