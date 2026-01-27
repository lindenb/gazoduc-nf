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

process JVARKIT_SVLEN {
	label "process_single"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	afterScript "rm -rf TMP"
	tag "${meta.id}"
	input:
		tuple val(meta  ),path(vcf)
	output:
		tuple val(meta ),path("*.vcf.gz"),emit:vcf
		path("versions.yml"),emit:versions
	script:
	    def jvm = " -XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"
        def args1 = task.ext.args1?:""
        def args2 = task.ext.args2?:""
        def args3 = task.ext.args3?:""
        def prefix =  task.ext.prefix?:"${meta.id}.svlen"
	"""
	mkdir -p TMP

    echo 'if(variant.hasAttribute("SVLEN")) return true; return new VariantContextBuilder(variant).attribute("SVLEN",variant.getLengthOnReference()).make();' > TMP/jeter.code

	bcftools view  ${args1} "${vcf}" |\\
    awk 'BEGIN {
            found=0;
            }
        /^##INFO=<ID=SVLEN,/ {
            found = 1;
            print;
            next;
            }
        /^#CHROM/ {
            if(found==0) {
                printf("##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\\"Difference in length between REF and ALT alleles\\">\\n");
                }
            print;
            next;    
            }
            {
            print;
            }' |\
	jvarkit ${jvm} vcffilterjdk \\
        ${args2} \\
        --script  TMP/jeter.code | \\
	bcftools view ${args3} -O z -o TMP/jeter.vcf.gz

mv -v TMP/jeter.vcf.gz "${prefix}.vcf.gz"

cat << END_VERSIONS > versions.yml
"${task.process}":
	jvarkit: todo
END_VERSIONS
	"""


stub:
	def prefix = "${meta.prefix}.svlen"
"""
touch versions.yml ${prefix}.vcf.gz 
"""
}
