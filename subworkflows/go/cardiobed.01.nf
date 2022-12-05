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

include {getVersionCmd;jvarkit;isHg19;isHg38;isBlank;hasFeature;moduleLoad;getGnomadGenomePath;getGnomadExomePath} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {GO_TO_BED_01} from './go2bed.01.nf'

/** CALL GO_TO_BED_01 for cardio-related terms */
workflow CARDIOBED_01 {
	take:
		meta
		reference
		gff3
	main:
		version_ch = Channel.empty()

		go_ch = GO_TO_BED_01(meta.plus(
			"goTerms" : "GO:0034765,GO:0043269,GO:0005216,GO:0006811,GO:0034220,GO:0008016",
			"excludeGoTerms" : "GO:0045202,GO:0099536,GO:0060078,GO:0003014" 
			), reference, gff3)
                version_ch = version_ch.mix(go_ch.version)

		version_ch = MERGE_VERSION(meta, "cardiobed", "Bed for Cardiac Genes", version_ch.collect())
	emit:
		version = version_ch
		obo = go_ch.obo
		goa = go_ch.goa
		bed = go_ch.bed
		index = go_ch.index
		terms = go_ch.terms
	}
