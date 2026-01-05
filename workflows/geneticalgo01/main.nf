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
nextflow.enable.dsl=2



include {runOnComplete;moduleLoad;dumpParams} from '../../modules/utils/functions.nf'
include {GENETICALGO as STEP1} from '../../subworkflows/geneticalgo01/geneticalgo01.nf' addParams([step_id:"step1",step_name:"Iteration 1"])
include {GENETICALGO as STEP2} from '../../subworkflows/geneticalgo01/geneticalgo01.nf' addParams([step_id:"step2",step_name:"Iteration 2"])
include {GENETICALGO as STEP3} from '../../subworkflows/geneticalgo01/geneticalgo01.nf' addParams([step_id:"step3",step_name:"Iteration 3"])
include {MULTIQC} from '../../subworkflows/multiqc/multiqc.nf'

if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}



workflow {
        RUN_GENALGO(
		params.genomeId,
		file(params.vcf),
		file(params.cases),
		file(params.controls)
		)
        }

runOnComplete(workflow)


workflow RUN_GENALGO {
	take:
		genomeId
		vcf
		cases
		controls
	main:
		tsv_ch = Channel.empty()
		mqc_ch = Channel.empty()
	        ch1 = STEP1(genomeId,vcf,cases,controls,Channel.empty())
		tsv_ch = tsv_ch.mix(ch1.tsv)
		mqc_ch = mqc_ch.mix(ch1.mqc)

	        ch2 = STEP2(genomeId,vcf,cases,controls, ch1.exclude_bed)
		tsv_ch = tsv_ch.mix(ch2.tsv)
		mqc_ch = mqc_ch.mix(ch2.mqc)

	        ch3 = STEP3(genomeId,vcf,cases,controls, ch2.exclude_bed)
		tsv_ch = tsv_ch.mix(ch3.tsv)
		mqc_ch = mqc_ch.mix(ch3.mqc)

		plot_ch=PLOT(tsv_ch.map{T->T[0].selname+" "+T[0].step_name+"\t"+T[1]}.collect())
		mqc_ch = mqc_ch.mix(plot_ch.mqc)
	
		mqc_ch = MULTIQC(mqc_ch.flatten().collect())

        }

process PLOT {
executor "local"
tag "N=${L.size()}"
input:
	val(L)
output:
	path("plot_mqc.yml"),emit:mqc
script:
"""
mkdir -p TMP
export LC_ALL=C


cat << EOF > TMP/jeter.tsv
${L.join("\n")}
EOF

cat << EOF > TMP/jeter.yml
id: "my_pca_section4"
parent_id: custom_section
parent_name: "Digest"
parent_description: "All steps are summarized here"
section_name: "PValue per generation"
description: "Plot each pvalue per generation. WARNING FOR NOW THERE IS A PROBLEM WITH THE X AXIS I NEED TO FIX IT"
plot_type: "linegraph"
pconfig:
      id: "example_coverage_lineplot"
      title: "PValue=f(Generations)"
      xlab: "Generations"
      ylab: "-log10(pvalue)"
data:
EOF

cat TMP/jeter.tsv | while IFS="\t" read N F
do
	awk -F '\t' -vN="\${N}" '(NR==1) {for(i=0;i<=NF;i++) {PREV=-1;if(\$i=="GENERATION") {X=i;} else if(\$i=="PVALUE") Y=i;}next;} {if(NR==2) {printf("  \\"%s\\":\\n",N); } if(int(\$X)==PREV) next;PREV=int(\$X); printf("    \\"%s\\": %f\\n",\$X,-(log(\$Y)/log(10)));}' "\$F" >> TMP/jeter.yml
done

mv TMP/jeter.yml plot_mqc.yml

"""
}
