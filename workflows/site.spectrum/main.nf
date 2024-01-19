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
nextflow.enable.dsl=2

/*

Allele frequency spectrum

> In population genetics, the allele frequency spectrum, sometimes called the site frequency spectrum, is the distribution of the allele frequencies of a given set of loci (often SNPs) in a population or sample

*/

include {moduleLoad;dumpParams;runOnComplete} from '../../modules/utils/functions.nf'
include {MULTIQC} from '../../subworkflows/multiqc/multiqc.nf'

if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}

workflow {
	ch1 = AF_SPECTRUM(params.genomeId, file(params.vcf), file(params.sample2population))
	html = VERSION_TO_HTML(ch1.version)
	}

workflow AF_SPECTRUM {
	take:
		genomeId
		vcf
		sample2population
	main:
		pop_ch = sample2population.splitCsv(sep:'\t',header:false).map{T->T[1]}.unique()
		by_ch = BY_POP(genomeId,vcf,sample2population,pop_ch)

		pairs_ch = by_ch.output.
			combine(by_ch.output).
			filter{T->T[0].compareTo(T[2])<0}
		
		plot_ch = PLOT(vcf,pairs_ch)

		mqc_ch = MULTIQC(plot_ch.output.collect())
}


process BY_POP {
tag "${pop}"
input:
	val(genomeId)
	path(vcf)
	path(sample2population)
	val(pop)
output:
	tuple val(pop),path("${pop}.count"),emit:output
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
"""
hostname 1>&2
${moduleLoad("bcftools")}
set -o pipefail
mkdir -p TMP

awk -F '\t' '(\$2=="${pop}") {print \$1;}' '${sample2population}' | sort -T TMP | uniq > TMP/jeter.samples.txt
test -s TMP/jeter.samples.txt

bcftools ${vcf.name.endsWith(".list")?"merge -a ":"view" -O u "${vcf}" |\
	bcftools view -m2 -M2 --samples-file  TMP/jeter.samples.txt --trim-alt-alleles -O u |\
	bcftools norm -i 'AC>0' --fasta-ref "${reference}" --multiallelics -both -Ou |\
	bcftools query -i 'AC>0' -f '%INFO/AC\\n' |\
	sort -T TMP -t '\t' -k1,1 |\
	uniq -c |\
	awk '{printf("%s\t%s\\n",\$2,\$1);}' |\
	sort -T TMP -t '\t' -k1,1  > TMP/jeter.count

mv TMP/jeter.count "${pop}.count"

"""
}

process PLOT {
tag '${pop1} vs ${pop2}"
input:
	path(vcf)
	tuple val(pop1),path(count1),val(pop2),path(count2)
output:
	path("${pop1}_${pop2}.yaml"),emit:output
script:
"""

cat << EOF > jeter.yaml
id: "my_sfs_${pop1}_${pop2}"
section_name: "SFS ${pop1} ${pop2}"
description: "Site Frequencey Spectrum for population ${pop1} and ${pop2} in vcf ${vcf}."
plot_type: "scatter"
pconfig:
  id: "sfs_${pop1}_${pop2}"
  title: "Site frequency spectrum"
  xlab: "${pop1}"
  ylab: "${pop2}"
data:
EOF

join -t '\t' -1 1 -2 1 -a 0 -a 1 -e 0 -o '1.1,1.2,2.2' "${count1}" "${count2}" |\
	awk '{printf("  AC_%s: { x: %d, y: %d }\\n,",\$1,\$2,\$3);}' >> jeter.yaml


mv jeter.yaml "${pop1}_${pop2}.yaml"
"""
}



runOnComplete(workflow);

