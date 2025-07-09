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

include {VCF_TO_BED} from '../vcf2bed'
include {BCFTOOLS_CONCAT_ALL} from '../bcftools/concat.all'


Map groupByFunction(hash,meta) {
	return [id:"truvari"];
}

workflow TRUVARI {
     take:
		meta
		fasta
		fai
		dict
        vcfs /* [meta,vcf,idx] */
     main:
		APPLY_TRUVARI(
			fasta,
			fai,
			dict,
			vcfs
				.map{[groupByFunction(it,meta), [it[1],it[2]] ]}
				.groupTuple()
				.map{[it[0],it[1].flatten()]}
		)

	emit:
		vcf = APPLY_TRUVARI.out.vcf
	}

process APPLY_TRUVARI {
    label "process_short"
	tag "${meta.id}"
    afterScript "rm -rf TMP"
    conda "${moduleDir}/../../conda/truvari.01.yml"
    input:
		tuple val(meta1),path(fasta)
		tuple val(meta2),path(fai)
		tuple val(meta3),path(dict)
		tuple val(meta ),path("VCFS/*")
    output:
		tuple val(meta),path("*.bcf"),path("*.bcf.csi"),emit: vcf
		path("versions.yml"),emit:versions
    script:
		def args1 = task.ext.args1?:""
		def args2 = task.ext.args2?:""
		def args3 = task.ext.args3?:""

		def prefix= task.ext.prefix?:meta.id
    """
	hostname 1>&2
	mkdir -p TMP

	find VCFS -type l \\( -name "*.vcf.gz" -o -name "*.bcf" \\) |\\
		while read F
		do
			bctools query -l "\${F}" | head -n 1 | tr "\\n" "\t" >> TMP/jeter.txt
			echo "\${F}" >> TMP/jeter.txt
		done
	
	sort -T TMP -t '\t' -k1,1 TMP/jeter.txt | cut -f 2 > TMP/jeter.list

	# bug FORMAT pour dragen
	bcftools merge --threads ${task.cpus}  ${args1} --force-samples --filter-logic '+'  --file-list TMP/jeter.list  -m none -O u -o TMP/jeter2.bcf
	mv TMP/jeter2.bcf TMP/jeter.bcf

	# optional filter ?
	bcftools view --threads ${task.cpus} ${args2}  -O u -o TMP/jeter2.bcf TMP/jeter.bcf
	mv TMP/jeter2.bcf TMP/jeter.bcf

	# bug DRAGEN
	bcftools annotate --threads ${task.cpus} --force -x 'FORMAT/SR' -O z -o TMP/merged.vcf.gz TMP/jeter.bcf
	bcftools index --threads ${task.cpus}  --tbi TMP/merged.vcf.gz

	# invoke truvari
	truvari collapse ${args3} --reference "${fasta}" -i "TMP/merged.vcf.gz" -c TMP/collapsed.vcf.gz |\\
		bcftools view -O u -o TMP/jeter.bcf
	
	bcftools +fill-tags --threads ${task.cpus}  -O u -o TMP/jeter2.bcf TMP/jeter.bcf -- -t AN,AC,AF 
	mv TMP/jeter2.bcf TMP/jeter.bcf
		
	bcftools sort -T TMP/sort -O b -o TMP/${prefix}.bcf TMP/jeter.bcf

	bcftools index --threads ${task.cpus} -f TMP/${prefix}.bcf
	mv TMP/${prefix}.bcf ./
	mv TMP/${prefix}.bcf.csi ./


cat << END_VERSIONS > versions.yml
"${task.process}":
	truvari: todo
END_VERSIONS
    """
   }
