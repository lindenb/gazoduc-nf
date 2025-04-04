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


include {TRUVARI_01} from '../../truvari'

workflow MANTA_SINGLE_SV01 {
	take:
		samplesheet
		bed
	main:
		version_ch = Channel.empty()

		sn_bam_bai_ch = samplesheet.
			splitCsv(header:true,sep:'\t').
			map{T->{
				if(!T.containsKey("sample")) throw new IllegalArgumentException("sample missing in "+T);
				if(!T.containsKey("bam")) throw new IllegalArgumentException("bam missing in "+T);
				if(!T.containsKey("bai")) throw new IllegalArgumentException("bai missing in "+T);
				if(!T.containsKey("fasta") || T.fasta.equals(".") || T.fasta.isEmpty()) {
					T = T.plus(fasta:params.fasta, fai:params.fai)
					}
				return T;
				}}.map{[it.sample,it.bam,it.bai,it.fasta,it.fai]}

		manta_ch = APPLY_MANTA_SINGLE(sn_bam_bai_ch,bed)	
		
	
		
		if(params.with_merge_manta_vcf==false) {
			merge_vcf = Channel.empty()
			merge_vcf_index = Channel.empty()
			}
		else
			{
			/* TODO
			truvari_ch = TRUVARI_01([:] , genomeId, file_list_ch.output)
			to_zip = to_zip.mix(truvari_ch.version)

			merge_vcf = truvari_ch.vcf
			merge_vcf_index = truvari_ch.index

			*/
			}


		
		
	emit:
		merge_vcf = merge_vcf
		merge_vcf_index = merge_vcf_index
		//manta_files = file_list_ch.output
	}


process APPLY_MANTA_SINGLE {
    tag "${sample} ${bam.name}"
    label "process_quick"
    conda "${moduleDir}/../../../conda/manta.yml"
    afterScript "rm -rf TMP"
    input:
	tuple val(sample),path(bam),path(bai),path(fasta),path(fai)
	path(bed)
    output:
    	tuple val(sample),path("*.{vcf.gz,vcf.gz.tbi}"),emit:output
    script:
	"""
	set -x
	hostname 1>&2
	mkdir -p TMP

	if test -f "${bed}"
	then
		${bed.name.endsWith("gz")?"gunzip -c":"cat"} "${bed}" |\\
			cut -f1,2,3 |\\
			jvarkit bedrenamechr -R "${fasta}" --column 1 --convert SKIP |\\
			awk -F '\t' '(\$1 ~ /^${params.regex_contig}\$/ )' |\\
			LC_ALL=C sort -t '\t' -T TMP -k1,1 -k2,2n |\\
			bedtools merge > TMP/intervals.bed
	else
		awk -F '\t' '(\$1 ~ /^${params.regex_contig}\$/ ) {printf("%s\t0\t%s\\n",\$1,\$2);}' "${fai}" |\\
			LC_ALL=C sort -t '\t' -T TMP -k1,1 -k2,2n > TMP/intervals.bed
		
	fi
	
	test -s TMP/intervals.bed
	bgzip TMP/intervals.bed
	tabix -f -p bed TMP/intervals.bed.gz
	
	configManta.py  --bam "${bam}" --referenceFasta "${fasta}" \\
		--callRegions TMP/intervals.bed.gz \\
		--runDir "TMP"

	./TMP/runWorkflow.py --quiet -m local -j '${task.cpus}'
	
	rm -rf ./TMP/workspace

	
	# convert BND TO INVERSIONS (added 20230115 but not tested)
	DIPLOID=`find ./TMP -type f -name "diploidSV.vcf.gz"` 1>&2
	which configManta.py 1>&2
	test ! -z "\${DIPLOID}"
	\$(ls \$( dirname \$(which configManta.py) )/../share/manta*/libexec/convertInversion.py)  `which samtools` "${fasta}" "\${DIPLOID}" | bcftools sort -T TMP -O z -o TMP/jeter.vcf.gz

	bcftools index -t TMP/jeter.vcf.gz

	mv -v TMP/jeter.vcf.gz "\${DIPLOID}"
	mv -v TMP/jeter.vcf.gz.tbi "\${DIPLOID}.tbi"

	# not the same dictionary ? rename the vcf
	cmp "${fai}" "${params.fai}" ||  find ./TMP -type f -name "*.vcf.gz" | while read F
	do	
		TODO change dictionary
	done
	

	# change name to sample
	find ./TMP -type f -name "*.vcf.gz" \
		-printf 'mv -v %p  ${sample}.%f\\n mv -v %p.tbi ${sample}.%f.tbi\\n' |bash 
	
	ls *.vcf.gz
	"""
	}
