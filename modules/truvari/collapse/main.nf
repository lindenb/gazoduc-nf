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

/*


usage: collapse [-h] -i INPUT [-o OUTPUT] [-c REMOVED_OUTPUT] [-f REFERENCE] [-w] [-k {first,maxqual,common}] [--bed BED] [--gt {off,all,het}]
                [--intra] [--median-info] [-F] [--debug] [-r REFDIST] [-p PCTSEQ] [-P PCTSIZE] [-O PCTOVL] [-t] [-n] [-m MAX_RESOLVE] [-D] [-d]
                [-B BNDDIST] [--hap] [--chain] [--no-consolidate] [-s SIZEMIN] [-S SIZEMAX] [--passonly]



*/

process TRUVARI_COLLAPSE {
    label "process_short"
	tag "${meta.id}"
    afterScript "rm -rf TMP"
    conda "${moduleDir}/../../../conda/truvari.01.yml"
    input:
		tuple val(meta1),path(fasta)
		tuple val(meta2),path(fai)
		tuple val(meta4),path(bed) //see https://github.com/ACEnglish/truvari/wiki/Collapse-on-Regions-with-a-High%E2%80%90Density-of-SVs
		tuple val(meta ),path(vcf),path(tbi)//generated with bcftools merge
    output:
		tuple val(meta),path("*.merged.vcf.gz"),path("*.merged.vcf.gz.tbi"),emit: vcf
		tuple val(meta),path("*.collapsed.vcf.gz"),emit: collapsed
		path("versions.yml"),emit:versions
    script:
		def args1 = task.ext.args1?:""
		def prefix= task.ext.prefix?:meta.id
		def merge_mode = task.ext.merge_mod?:"id"
    """
	hostname 1>&2
	mkdir -p TMP
	

	# invoke truvari
	truvari collapse \\
			${args1} \\
			${bed?"--bed \"${bed}\"":""} \\
			--reference "${fasta}" \\
			-i "${vcf}" \\
			-c TMP/collapsed.vcf.gz |\\
		bcftools view -O u -o TMP/jeter.bcf
	
	bcftools +fill-tags \\
		--threads ${task.cpus} \\
		-O u \\
		-o TMP/jeter2.bcf TMP/jeter.bcf -- -t  AN,AC,AF,AC_Hom,AC_Het,AC_Hemi,NS
	mv TMP/jeter2.bcf TMP/jeter.bcf
	
	bcftools sort -T TMP/sort -O z -o TMP/${prefix}.merged.vcf.gz TMP/jeter.bcf

	bcftools index -t --threads ${task.cpus} -f TMP/${prefix}.merged.vcf.gz

	mv TMP/${prefix}.merged.vcf.gz ./
	mv TMP/${prefix}.merged.vcf.gz.tbi ./

	mv TMP/collapsed.vcf.gz ${prefix}.collapsed.vcf.gz


cat << END_VERSIONS > versions.yml
"${task.process}":
	truvari: \$(truvari version)
END_VERSIONS
    """

stub:
def prefix= task.ext.prefix?:meta.id
"""
touch versions.yml ${prefix}.merged.vcf.gz ${prefix}.merged.vcf.gz.tbi ${prefix}.collapsed.vcf.gz
"""
}
