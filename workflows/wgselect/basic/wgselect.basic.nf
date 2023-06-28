/*

Copyright (c) 2023 Pierre Lindenbaum

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

def gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults()



gazoduc.
    make("vcf","NO_FILE").
    description("indexed vcf").
    required().
    put()

gazoduc.
    make("pedigree","NO_FILE").
    description("pedigree").
    put()

gazoduc.
    make("bed","NO_FILE").
    description("pedigree").
    put()

include {runOnComplete;moduleLoad;getKeyValue;hasFeature} from '../../../modules/utils/functions.nf'
include {WGSELECT_02} from '../../../subworkflows/wgselect/wgselect.02.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'

if(params.help) {
        gazoduc.usage().name("wgselect").description("Cleanup variants from a VCF").print()
        exit 0
	}
else
    	{
	gazoduc.validate()
        }


workflow {
		ch1 = WGSELECT_02(params.genomes[params.genomeId], file(params.vcf), file(params.pedigree), file(params.bed))
		html = VERSION_TO_HTML(ch1.version)
		PUBLISH(params, ch1.contig_vcfs, ch1.variants_list , ch1.version, html.html)
		}

runOnComplete(workflow)

process PUBLISH {
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
executor "local"
input:
	val(meta)
	path(vcfs)
	path(variants)
	path(xml)
	path(html)
output:

when:
        !params.getOrDefault("publishDir","").trim().isEmpty()
script:
	prefix = meta.prefix?:""
"""
mkdir "${prefix}archive"

for F in "${vcfs}" "${variants}" "${xml}" "${html}"
do
	ln -s "\${F}" "./${prefix}archive/" 
done

zip -r "${prefix}archive.zip" "${prefix}archive"
"""
}
