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
process CALL_DELLY {
    tag "${meta.id}"
    label "process_short"
    afterScript "rm -rf TMP"
    conda "${moduleDir}/../../../conda/delly.yml"
    input:
		tuple val(meta1),path(fasta)
		tuple val(meta2),path(fai)
		tuple val(meta3),path(optional_exclude_bed)
		tuple val(meta4),path(optional_genotype_vcf),path(optional_genotype_vcf_tbi)
		tuple val(meta),path(bam),path(bai)
    output:
         	tuple val(meta),path("*.bcf"),path("*.bcf.csi"),emit:vcf
		path("versions.yml"),emit:versions
    script:
		def args1  = task.ext.args1?:""
		def prefix = task.ext.prefix?:meta.id
	"""
	hostname 1>&2
	mkdir -p TMP
	export TMPDIR=\${PWD}/TMP
	export PATH=\${PWD}:\${PATH}

	delly call \\
		${args1} \\
		${optional_exclude_bed?"--exclude \"${optional_exclude_bed}\"":""} \\
		${optional_genotype_vcf?"--vcffile \"${optional_genotype_vcf}\"":""} \\
		--outfile "TMP/${prefix}.bcf" \\
		--genome "${fasta}" \\
		"${bam}" 1>&2

	mv -v TMP/${prefix}.* ./
	
	touch versions.yml
	"""

	stub:
		def prefix = task.ext.prefix?:meta.id
	"""
	touch ${prefix}.bcf ${prefix}.bcf.csi versions.yml
	"""
    }
