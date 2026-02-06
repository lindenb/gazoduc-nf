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
include {INDEXCOV_REBUILD_BAI                } from '../../../modules/indexcov/rebuild.bai'
include {APPLY_INDEXCOV                      } from '../../../modules/indexcov/apply'
include {MERGE_BEDS                          } from '../../../modules/indexcov/merge.beds'
include {JVARKIT_INDEXCOV2VCF                } from '../../../modules/jvarkit/indexcov2vcf'
include {INDEXCOV_TO_SVG                     } from '../../../subworkflows/indexcov/circular'
include {SAMTOOLS_COVERAGE                   } from '../../../modules/samtools/coverage'
include {parseBoolean                        } from '../../../modules/utils/functions.nf'

workflow INDEXCOV {
     take:
	 	metadata
        fasta
		fai
		dict
		scatter_N_bed
		pedigree
		bams /* 		meta,bam,(bai)*/

     main:
        versions = Channel.empty()
		multiqc = Channel.empty()

		ch1 = bams.branch{
			rebuild: it.size()==2 || it[1].name.endsWith(".cram") || (it[0].force_rebuild && it[0].force_rebuild==true)
			no_need_to_rebuild: true
			}

		INDEXCOV_REBUILD_BAI(
			fasta,
			fai,
			ch1.rebuild.map{[it[0],it[1]]} //remove bai, if any
			)
		versions = versions.mix(INDEXCOV_REBUILD_BAI.out.versions)

		def collate_size = ((metadata.batch_size?:1000) as int)

		bams2_ch = INDEXCOV_REBUILD_BAI.out.bam.mix(ch1.no_need_to_rebuild)
			.map{[it[1],it[2]]} //bam,bai
			.collect(flat:false)
			.flatMap{
				def L1=[];
				def L2=[];
				def L3 = it.sort({a,b->a[0].name.compareTo(b[0].name)});
				int i=0;
				while(i < L3.size()) {
					if(L2.size()==collate_size*2 /* because bam and bai */) {
						L1.add(L2);
						L2=[]
						}
					L2.add(L3[i][0]);//bam
					L2.add(L3[i][1]);//bai
					i++;
					}
				if(L2.size()>0) L1.add(L2);
				return L1;
				}
			.map{T->[[id:"batch."+T.collect{V->V.toString()}.join(" ").md5().substring(0,7)], T]}
		
    	APPLY_INDEXCOV(fasta,fai,bams2_ch)
	versions = versions.mix( APPLY_INDEXCOV.out.versions)

	MERGE_BEDS(
		APPLY_INDEXCOV.out.bed.
			map{T->T[1]}.
			collect().
			map{[metadata,it]}
		)
	versions = versions.mix(MERGE_BEDS.out.versions)

	JVARKIT_INDEXCOV2VCF(
		fasta,fai,dict,
		pedigree,
		MERGE_BEDS.out.bed.map{meta,bed,_tbi->[meta,bed]}
		)
	versions = versions.mix(JVARKIT_INDEXCOV2VCF.out.versions)
	
	INDEXCOV_TO_SVG(
		metadata,
		fasta,
		fai,
		dict,
		scatter_N_bed,
		MERGE_BEDS.out.bed.map{meta,bed,_tbi->[meta,bed]}
		)
	versions = versions.mix(MERGE_BEDS.out.versions)

	if(metadata.with_singletons==null || parseBoolean(metadata.with_singletons)==true) {
	/* plot singletons */
	EXTRACT_SINGLETONS(
		MERGE_BEDS.out.bed
			.map{meta,bed,tbi->[meta,bed,bed]} //duplicate bed
			.splitCsv(header:false,sep:'\t',limit:1)
			.flatMap{meta,header,bed->{
				def min_n_samples = 20;
				def L=[];
				if(header.size()-3 < min_n_samples) {
					return L;
					}
				def i;
				for(i=3;i< header.size();i++) {
					L.add([
						[id:header[i],column0:i],
						bed
						])
					}
				return L;
				}}
			)
	versions = versions.mix(EXTRACT_SINGLETONS.out.versions)

	COLLECT_SINGLETONS(
		fai,
		EXTRACT_SINGLETONS.out.bed.map{meta,bed->bed}.collect().map{f->[[id:"singletons"],f.sort()]}
		)
	versions = versions.mix(COLLECT_SINGLETONS.out.versions)

	if(params.with_coverage==null || parseBoolean(params.with_coverage)==true) {
		SAMTOOLS_COVERAGE(
			fasta,
			fai,
			COLLECT_SINGLETONS.out.slop_bed
				.map{_meta,bed->bed}
				.splitCsv(sep:'\t',header:false)
				.map{row->[
						id : "${row[0]}_${row[1]}_${row[2]}",
						title : "${row[4]}",
						interval: "${row[0]}:${row[1]}-${row[2]}"
						]}
				.combine(
					ch1.rebuild
						.map{_meta,bam,bai->[bam,bai]}
						.flatMap()
						.collect()
						.map{f->[f.sort()]}
					)
			)
			
		versions = versions.mix(SAMTOOLS_COVERAGE.out.versions)
		}
	}

	

	emit:
		bed  = MERGE_BEDS.out.bed
		zip = APPLY_INDEXCOV.out.zip
		vcf  = JVARKIT_INDEXCOV2VCF.out.vcf
		versions
		multiqc
	}

process EXTRACT_SINGLETONS {
label "process_single"
tag "${meta.id}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta),path(bed)
output:
	tuple val(meta),path("*.bed"),optional:true,emit:bed
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:"${meta.id}"
	def treshold = (task.ext.treshold?:0.1)
	def column0 = (meta.column0?:-1) as int
	def column1 = column0 +1
	def count_other = meta.count_other?:0
	def args1 = task.ext.args1?:" -d 17000"
"""
mkdir -p TMP

cat << '__EOF__' > TMP/jeter.awk
(NR==1) {
	next;
	}
	{
	N1=0;
	N2=0;
	for(i=4;i<=NF;i++) {
		v = (\$i * 1.0);
		is_sv =  ( v <= (0.5+${treshold}) || v >= (1.5-${treshold}) ?1:0);
		if(i==${column1}) {
			if(is_sv==0) next;
			N1=1;
			}
		else {
			if(is_sv==1) {
				N2++;
				}
			}
		}
	if(N1==1 && N2<=${count_other}) {
		printf("%s\t%s\t%d\\n",\$1,\$2,int(\$3)+1);
		}
	}
__EOF__

gunzip -c "${bed}" |\\
	awk -F '\t' -f TMP/jeter.awk |\\
	sort -T TMP -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n |\\
	bedtools merge ${args1} |\\
	sed 's|\$|\t${meta.id}|' > TMP/jeter.bed

if test -s TMP/jeter.bed
then
	mv -v  TMP/jeter.bed ${prefix}.bed
fi

cat <<  __EOF__ > versions.yml
${task.process}:
	awk: todo
__EOF__
"""
stub:
	def prefix=task.ext.prefix?:"${meta.id}"
"""
touch versions.yml ${prefix}.bed
"""
}


process COLLECT_SINGLETONS {
label "process_single"
tag "${meta.id}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fai)
	tuple val(meta ),path(beds)
output:
	tuple val(meta),path("*.slop.bed"),emit:slop_bed
	tuple val(meta),path("*.concat.bed"),emit:bed
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:"${meta.id}"
"""
mkdir -p TMP

cut -f1,2 "${fai}" |\\
    sort  -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter.genome


cat ${beds} |\\
	sort -T TMP -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n > TMP/jeter1.bed


awk -F '\t'  '{B=int(\$2);E=int(\$3);L=(E-B);printf("%s\t%d\t%d\t%s\t%s:%s-%s LEN=%d archetype:%s\\n",\$1,B,E,\$4,\$1,\$2,\$3,L,\$4);}' TMP/jeter1.bed |\\
	bedtools slop -i - -g TMP/jeter.genome -b 2.0 -pct |\\
	sort -T TMP -S ${task.memory.kilo} -t '\t' -k1,1 -k2,2n > TMP/jeter2.bed

mv TMP/jeter1.bed ${prefix}.concat.bed
mv TMP/jeter2.bed ${prefix}.slop.bed

cat <<  __EOF__ > versions.yml
${task.process}:
	awk: todo
__EOF__
"""
stub:
	def prefix=task.ext.prefix?:"${meta.id}"
"""
touch versions.yml ${prefix}.slop.bed ${prefix}.concat.bed
"""
}
