/*

Copyright (c) 2024 Pierre Lindenbaum

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

def gazoduc = gazoduc.Gazoduc.getInstance()



include {moduleLoad} from '../utils/functions.nf'

process RVTESTS01_VCF_01 {
tag "${file(vcf).name}"
afterScript "rm -rf TMP"
memory "10g"
input:
	val(meta)
	val(reference)
	val(vcf)
	val(pedigree)
output:
	path("assoc.list"),emit:output
	path("version.xml"),emit:version
script:
	def rvtest_arguments = params.rvtest.arguments?:""
"""
hostname 1>&2
${moduleLoad("rvtests bcftools jvarkit")}

mkdir -p TMP ASSOC

## https://github.com/zhanxw/rvtests/issues/80
cut -f 1 "${reference}.fai" > TMP/chroms.A.txt
sed 's/^chr//' TMP/chroms.A.txt > TMP/chroms.B.txt
paste TMP/chroms.A.txt TMP/chroms.B.txt > TMP/chroms.C.txt


bcftools annotate -O v --rename-chrs TMP/chroms.C.txt "${vcf}" |\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/vcfgenesplitter.jar \
		-E 'ANN/GeneId ANN/FeatureId VEP/GeneId VEP/Ensp VEP/Feature' \
		--manifest TMP/jeter.mf \
		-o \${PWD}/TMP

i=1
	awk -F '\t' '/^[^#]/ {S=\$4;gsub(/[\\/]+/,"_",S);G=\$5;gsub(/[\\/_\\-]+/,"_",G);K=\$6;gsub(/[\\/]+/,"_",K);printf("%s_%s_%s\t%s:%d-%d\t%s\\n",K,G,S,\$1,\$2,\$3,\$7);}' TMP/jeter.mf |\
	while read K RANGE V
	do
		bcftools index -t "\${V}"

		# build setFile
		echo "\$K\t\${RANGE}" > TMP/variants.setfile

		rvtest  --noweb \
        		--inVcf "\${V}" \
			--setFile TMP/variants.setfile \
	        	--pheno "${pedigree}" \
		        --out "ASSOC/part.\${i}" \
			${rvtest_arguments} 2> TMP/last.rvtest.log



		i=\$((i+1))
	done


find \${PWD}/ASSOC -type f -name "part.*assoc" > assoc.list

cat << EOF > version.xml
<properties id="${task.process}">
  <entry key="name">${task.process}</entry>
  <entry key="description">VCF is split by transcript / gene, and then rvtest is applied.</entry>
  <entry key="name"Pedigree">${pedigree}</entry>
  <entry key="rvtest.path">\$(which rvtest)</entry>
  <entry key="rvtest.version">\$(rvtest  --version 2> /dev/null  | head -n1)</entry>
</properties>
EOF
"""
stub:
"""
touch assoc.list
echo "<properties/>" > version.xml
"""
}
