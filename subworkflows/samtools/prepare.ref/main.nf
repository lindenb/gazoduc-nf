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
include {SAMTOOLS_DICT         } from '../../../modules/samtools/dict'
include {SAMTOOLS_FAIDX        } from '../../../modules/samtools/faidx'
include {FAI2BED               } from '../../../modules/samtools/fai2bed'
include {SCATTER_TO_BED        } from '../../../subworkflows/gatk/scatterintervals2bed'
include {BEDTOOLS_COMPLEMENT   } from '../../../modules/bedtools/complement'
include {k1_signature          } from '../../../modules/utils/k1.nf'

/**
 * from a FASTA Channel, create fai, dict, scatter (bed for ACGT bases only )
 */
workflow PREPARE_REFERENCE {
take:
	workflow_meta
	fasta //meta,fasta
main:
	versions = Channel.empty()
	multiqc = Channel.empty()

	SAMTOOLS_FAIDX(fasta)
	versions = versions.mix(SAMTOOLS_FAIDX.out.versions)
	
	FAI_TO_BUILD(SAMTOOLS_FAIDX.out.fai)
	versions = versions.mix(FAI_TO_BUILD.out.versions)

	build_ch  = FAI_TO_BUILD.out.csv
		.splitCsv(header:true,sep:',')
		.map{[it.id,it]}

	fasta = fasta
		.map{meta1,fasta->[meta1.id,meta1,fasta]}
		.join(build_ch)
		.map{id,meta1,fasta,meta2->[meta1.plus(meta2),fasta]}
		
	
	fai = SAMTOOLS_FAIDX.out.fai
		.map{meta1,fai->[meta1.id,meta1,fai]}
		.join(build_ch)
		.map{id,meta1,fai,meta2->[meta1.plus(meta2),fai]}

	
	SAMTOOLS_DICT(fasta)
	versions = versions.mix(SAMTOOLS_DICT.out.versions)
	
	FAI2BED(fai)
	versions = versions.mix(FAI2BED.out.versions)
	
	
	complement_bed =  Channel.of([[id:"nocomplement"],[]])
	if(workflow_meta.skip_scatter==null || workflow_meta.skip_scatter==false) {
		SCATTER_TO_BED(
			workflow_meta,
			fasta,
			SAMTOOLS_FAIDX.out.fai,
			SAMTOOLS_DICT.out.dict
			)
		versions = versions.mix(SCATTER_TO_BED.out.versions)
		scatter_bed = SCATTER_TO_BED.out.bed

		if(workflow_meta.skip_complement==null || workflow_meta.skip_complement==false) {
			ch1 = fai.map{meta,fai->[meta.id,meta,fai]}
				.join(SCATTER_TO_BED.out.bed.map{meta,bed->[meta.id,meta,bed]})
				.multiMap{id,meta1,fai,meta2,bed->
					fai: [meta1,fai]
					bed:  [meta2,bed]
				}
			BEDTOOLS_COMPLEMENT(
				ch1.fai,
				ch1.bed
				)
			complement_bed = BEDTOOLS_COMPLEMENT.out.bed
			versions = versions.mix(BEDTOOLS_COMPLEMENT.out.versions)
			}
		}
	else
		{
		scatter_bed = Channel.of([[id:"noscatterbed"],[]])
		}

emit:
	versions
	fasta
	fai
	multiqc
	dict = SAMTOOLS_DICT.out.dict
	bed = FAI2BED.out.bed
	scatter_bed
	complement_bed
}

process FAI_TO_BUILD {
label "process_single"
tag "${meta.id} ${fai.name}"
conda "${moduleDir}/../../../conda/bioinfo.02.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta),path(fai)
output:
	path("*.csv"),emit:csv
	path("versions.yml"),emit:versions
script:	
	def k1 = k1_signature()
    def prefix = task.ext.prefix?:"${meta.id}"
"""
mkdir -p TMP

echo "id"          > TMP/jeter1.txt
echo "${meta.id}"  > TMP/jeter2.txt

echo "ucsc_name" >> TMP/jeter1.txt

cat << EOF | sort -t '\t' -k1,1 > TMP/jeter.a
1:${k1.hg19}\thg19\t9606
1:${k1.hg38}\thg38\t9606
1:${k1.canFam3}\tcanFam3\t9612
1:${k1.canFam4}\tcanFam4\t9612
RF01:3302\trotavirus_rf\t10912
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' |\\
	sed 's/^chr//' |\\
	sort -t '\t' -k1,1 -T TMP > TMP/jeter.b

join -t '\t' -1 1 -2 1 -o 1.2 TMP/jeter.a TMP/jeter.b |\\
	awk -F '\t' 'BEGIN {B="undefined"} {B=\$1;} END {printf("%s\\n",B);}' >> TMP/jeter2.txt

echo "taxon_id" >> TMP/jeter1.txt
join -t '\t' -1 1 -2 1 -o 1.3 TMP/jeter.a TMP/jeter.b |\\
	awk -F '\t' 'BEGIN {B=""} {B=\$1;} END {printf("%s\\n",B);}' >> TMP/jeter2.txt

echo "chr_prefix" >> TMP/jeter1.txt
awk -F '\t' 'BEGIN {B="false"} /^chr/{B="true";} END {printf("%s\\n",B);}' '${fai}' >> TMP/jeter2.txt

paste -sd, TMP/jeter1.txt >  TMP/samplesheet.csv
paste -sd, TMP/jeter2.txt >> TMP/samplesheet.csv

mv TMP/samplesheet.csv "${prefix}.samplesheet.csv"
touch versions.yml
"""
}
