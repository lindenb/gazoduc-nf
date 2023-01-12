/*

Copyright (c) 2022 Pierre Lindenbaum

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


def gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults().putReference()

gazoduc.make("controls","").
	description("file containing the path to indexed BAM or CRAM files for controls").
	argName("file").
	put()

gazoduc.make("cases","").
	description("file containing the path to indexed BAM or CRAM files for cases").
	required().
	existingFile().
	put()

gazoduc.make("bnd",true).
	description("keep BND data").
	menu("delly").
	setBoolean().
	put()

gazoduc.make("cnv",true).
	description("run 'delly cnv'").
	menu("delly").
	setBoolean().
	put()


include {DELLY2_RESOURCES} from '../../subworkflows/delly2/delly2.resources.nf' 
include {DELLY2_SV} from '../../subworkflows/delly2/delly2.sv.nf' 
include {SAMTOOLS_CASES_CONTROLS_01} from '../../subworkflows/samtools/samtools.cases.controls.01.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'

if(params.help) {
	gazoduc.usage().
		name("delly2").
		description("Detects CNV/SV using delly2 ( https://github.com/dellytools/delly ) ").
		print();
	exit 0
	}
else	{
	gazoduc.validate()
	}


workflow {
	delly_ch = DELLY2_SV(params, params.reference, params.cases, params.controls)


	html = VERSION_TO_HTML(params,delly_ch.version)	

	to_publish = Channel.empty()
	to_publish = to_publish.
			mix(delly_ch.sv_vcf).
			mix(delly_ch.sv_vcf_index).
			mix(delly_ch.cnv_vcf).
			mix(delly_ch.cnv_vcf_index).
			mix(delly_ch.version).
			mix(html.html)

	PUBLISH(to_publish.collect())
	}

process PUBLISH {
tag "N=${files.size()}"
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	val(files)
output:
	path("*.bcf")
	path("*.bcf.csi")
	path("*.xml")
	path("*.html")
when:
	!params.getOrDefault("publishDir","").trim().isEmpty()
script:
	prefix = params.getOrDefault("prefix","")
"""
for F in ${files.join(" ")}
do
	ln -s "\${F}" ./
done
"""
}

workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}

