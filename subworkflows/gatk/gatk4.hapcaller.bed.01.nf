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
include {GATK4_HAPCALLER_DIRECT_01} from './gatk4.hapcaller.direct.01.nf'
include {BED_CLUSTER_01} from '../../modules/jvarkit/jvarkit.bedcluster.01.nf'
include {getKeyValue} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'


workflow GATK4_HAPCALLER_BED_01 {
	take:
		meta
		reference
		references
		bams
		bed
		dbsnp
		pedigree
	main:
		def max_bams_direct = (getKeyValue(meta,"max_bams_gatk_direct","20") as Integer);
		version_ch = Channel.empty()
		
		cluster_ch = BED_CLUSTER_01(meta.plus(
			["bed_cluster_method":getKeyValue(meta,"bed_cluster_method","--size 1mb")]),
			reference,
			bed
			)
		cluster_ch.output.view{"Clustering list of beds is : $it ."}
		version_ch = version_ch.mix(cluster_ch.version)
		
		
		gatk_strategy = Channel.fromPath(bams).splitText().count().combine(Channel.fromPath(bams)).branch{
			direct: (it[0] as int) <= max_bams_direct
			gvcfs:  (it[0] as int) > max_bams_direct
			}

		gatk_ch = GATK4_HAPCALLER_DIRECT_01(
					meta.plus(["nbams":max_bams_direct]),
					reference,
					references,
					gatk_strategy.direct.map{T->T[1]},
					cluster_ch.output,
					dbsnp,
					pedigree)
		version_ch = version_ch.mix(gatk_ch.version)

		version_ch = MERGE_VERSION(meta, "gatk4 bed", "call bams in a bed region", version_ch.collect())

	emit:
		vcf = gatk_ch.vcf
		index = gatk_ch.index
		version = version_ch.version
	}

