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

include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {getKeyValue;moduleLoad;getVersionCmd} from '../../modules/utils/functions.nf'
include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {BCFTOOLS_CONCAT_01} from '../bcftools/bcftools.concat.01.nf'
include {SURVIVOR_INSTALL} from '../../modules/survivor/survivor.install.nf'
include {SURVIVOR_MERGE} from '../../modules/survivor/survivor.merge.nf'

def gazoduc = gazoduc.Gazoduc.getInstance(params)

/**

merge VCFs using survivor

*/
workflow SURVIVOR_01 {
     take:
        meta /* meta */
        reference /* indexed fasta reference */
        vcfs /* file containing the path to VCF files . One per line */
     main:
        version_ch = Channel.empty()

	vcf2bed_ch = VCF_TO_BED(meta,vcfs)
	version_ch = version_ch.mix( vcf2bed_ch.version)

	exe_ch = SURVIVOR_INSTALL(meta)
	version_ch = version_ch.mix(exe_ch.version)


	perctg_ch = SURVIVOR_MERGE(meta, exe_ch.executable, vcf2bed_ch.chromosomes.splitText().map{it.trim()},vcfs)
	version_ch = version_ch.mix(perctg_ch.version)


	x3_ch = COLLECT_TO_FILE_01([:], perctg_ch.vcf.collect())
	version_ch = version_ch.mix(x3_ch.version)


	concat_ch = BCFTOOLS_CONCAT_01([:], x3_ch.output)
	version_ch = version_ch.mix(concat_ch.version)
		
		
	version_ch = MERGE_VERSION(meta, "Survivor", "Survivor Merge SV Vcfs",version_ch.collect())

	emit:
		vcf = concat_ch.vcf
		index = concat_ch.index
		version= version_ch
	}
