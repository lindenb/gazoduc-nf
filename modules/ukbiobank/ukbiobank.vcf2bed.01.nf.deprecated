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

include {moduleLoad;parseBoolean} from '../utils/functions.nf'

process UKBIOBANK_VCF2BED_01 {
	executor "local"
	input:
		val(meta)
	output:
		path("ukbiobank.vcf.bed"),emit:bed
		path("version.xml"),emit:version
	script:
		def base = "/LAB-DATA/BiRD/shares/ITX/u1087/ukbiobank_49823"
		def blocks = "${base}/exomes/OQFE_pVCF_23156/pvcf_blocks.txt"
	"""
	hostname 1>&2
	set -o pipefail
	
	awk -F '\t' '{C=\$2;if(C=="23")C="X";if(C=="24") C="Y";printf("chr%s\t%s\t%s\t${base}/exomes/OQFE_pVCF_23156/ukb23156_c%s_b%s_v1.vcf.gz\\n",C,\$4,\$5,C,\$3);}' "${blocks}" > ukbiobank.vcf.bed


	#####################################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Convert biobank vcf to BED</entry>
		<entry key="blocks">${blocks}</entry>
	</properties>
	EOF
	"""
	}
