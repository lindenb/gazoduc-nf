/*

Copyright (c) 2025 Pierre Lindenbaum

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



include {runOnComplete               } from '../../modules/utils/functions.nf'
include {JVARKIT_BAM_WITHOUT_BAI     } from '../../modules/jvarkit/bamwithoutbai'
include {BAM_QC                      } from '../../subworkflows/bamqc'
include {MULTIQC                     } from '../../modules/multiqc'
include {COMPILE_VERSIONS            } from '../../modules/versions'
include {BATIK_DOWNLOAD              } from '../../modules/batik/download'
include {BATIK_RASTERIZE             } from '../../modules/batik/rasterize'
include {GHOSTSCRIPT_MERGE           } from '../../modules/gs/merge'
runOnComplete(workflow)

workflow {
	def hash_ref= [
		id: file(params.fasta).baseName,
		name: file(params.fasta).baseName,
		ucsc_name: (params.ucsc_name?:"undefined"),
		ensembl_name: "GRCh38"
		]
	def fasta = [ hash_ref, file(params.fasta)]
	def fai   = [ hash_ref, file(params.fai)]
	def dict  = [ hash_ref, file(params.dict)]
	def gtf = [hash_ref,file(params.gtf),file(params.gtf+".tbi")]
	def bed  = [ hash_ref, file(params.bed)]
	if(!hash_ref.ucsc_name.equals("hg38")) throw new IllegalArgumentException("WANT hg38");
	

	version_ch = Channel.empty()
	multiqc = Channel.empty()
		//batik_ch = BATIK_DOWNLOAD_01(meta)
		//version_ch = version_ch.mix(batik_ch.version)


	ENCODE_METADATA(hash_ref)
	version_ch = version_ch.mix(ENCODE_METADATA.out.versions)

	ch1 = ENCODE_METADATA.out.tsv.map{it[1]}.splitCsv(header:true,sep:'\t')
		.filter{it.File_format.equals("bam")}
		.filter{it.Output_type.equals("alignments")}
		.filter{it.File_assembly.equals("GRCh38")}
		.filter{!it.File_accession.isEmpty()}
		.filter{!it.File_download_URL.isEmpty()}
		.map{[[id:it.File_accession],it.File_download_URL]}
		.filter{!it[1].matches("ENCFF402BMB|ENCFF337XBW|ENCFF008ILD|ENCFF241TBK|ENCFF722YLK|ENCFF818YGY|ENCFF849EUO|ENCFF979GXG") }//broken bam
	
	
	ch2 = Channel.fromPath(params.bed).splitCsv(header: false,sep:'\t',strip:true)
			.map{[
				contig:it[0],
				start:(it[1] as int) +1,
				end: it[2] as int,
				title: it.size() > 3 ? it[3] : "${it[0]}:${(it[1] as int)+1}-${it[2]}"
				]}
	
	FETCH_GENES( gtf, ch2 )
	version_ch = version_ch.mix(FETCH_GENES.out.versions)

	JVARKIT_BAM_WITHOUT_BAI(
		dict,
		ch1.combine(ch2).map{
		[
		it[0].plus(it[2]),
		it[1]
		]})
	version_ch = version_ch.mix(JVARKIT_BAM_WITHOUT_BAI.out.versions)


	tosvg_ch = JVARKIT_BAM_WITHOUT_BAI.out.bam
		.combine(FETCH_GENES.out.gtf)
		.filter{it[0].contig.equals(it[3].contig)}
		.filter{(it[0].start as int)===(it[3].start as int)}
		.filter{(it[0].end as int)===(it[3].end as int)}
		.map{[it[0],it[1],it[2],it[4],it[5]]}
		
	
	PLOT_SVG(tosvg_ch)
	
		
	BATIK_DOWNLOAD(hash_ref);

	BATIK_RASTERIZE(
		BATIK_DOWNLOAD.out.rasterizer,
		PLOT_SVG.out.svg
		.map{[it[0],(it[1] instanceof List?it[1]:[it[1]])]}
			.flatMap{
				def L=[];
				for(int i=0;i< it[1].size();i++) {
					L.add([it[0],it[1][i]]);
					}
				return L;
				}
		)


	GHOSTSCRIPT_MERGE(BATIK_RASTERIZE.out.output
		.map{[[id:it[0].contig+"_"+it[0].start+"_"+it[0].end],it[1]]}
		.groupTuple()
	)
	

	ch1 = JVARKIT_BAM_WITHOUT_BAI.out.bam.map{[
			[id:it[0].id],
			[it[1],it[2]]//bam,bai
			]}
			.groupTuple()
			.map{[it[0],it[1].flatten()]}
			.branch{v->
				need_merge: v[1].size()>2//more than one bam and one bai
				other: true
				}

	

	ch1.need_merge.map{throw new IllegalStateException("TODO need to merge");}
	
	all_bams = ch1.other
		.map{[it[0],it[1][0],it[1][1]]}
	
	GET_EXONS(gtf,bed)
	version_ch = version_ch.mix(GET_EXONS.out.versions)

	BAM_QC(
		hash_ref,
		fasta,
		fai,
		dict,
		GET_EXONS.out.bed.first(),
		all_bams
		)
	version_ch = version_ch.mix(BAM_QC.out.versions)
	multiqc = multiqc.mix(BAM_QC.out.multiqc)


	// group BAMs in 255 pools using md5sum.subtr(0,2)

	pools_ch = all_bams
		.map{[[id:it[0].id.md5().substring(0,2)],[it[1],it[2]]]}
		.groupTuple()
		.map{[it[0],it[1].flatten()]}
	
	SCAN_FOR_SPLICE_EVENTS(pools_ch)
	CAT_SPLICE_EVENTS(SCAN_FOR_SPLICE_EVENTS.out.bed.map{[[id:"junctions"],it[1]]}.groupTuple())


	if(params.cryptic_tsv!=null) {
		SCAN_FOR_CRYPTIC(
			[[id:"cryptic"],file(params.cryptic_tsv)],
			pools_ch
			)
		CAT_CRYPTIC(SCAN_FOR_CRYPTIC.out.tsv.map{[[id:"cryptic"],it[1]]}.groupTuple())
	}
	

 	COMPILE_VERSIONS(version_ch.collect().map{it.sort()})
    multiqc = multiqc.mix(COMPILE_VERSIONS.out.multiqc.map{[[id:"versions"],it]})
    // in case of problem multiqc_ch.filter{!(it instanceof List) || it.size()!=2}.view{"### FIX ME ${it} MULTIQC"}
    MULTIQC(multiqc.map{it[1]}.collect().map{[[id:"sashimi.encode"],it]})

	/*	
		bam_ch = DOWNLOAD_BAM_01(meta,gtf, all_bams.combine(all_intervals))
		version_ch = version_ch.mix(bam_ch.version)

		merge_pdf = MERGE_PDF(meta,bam_ch.pdf.groupTuple())
		version_ch = version_ch.mix(merge_pdf.version)

		junction_ch = COLLECT_ALL_JUNCTIONS(meta,bam_ch.junctions.collect());
		version_ch = version_ch.mix(junction_ch.version)

		version_ch = MERGE_VERSION(meta, "encodeRNA", "sashimi plot encode RNA", version_ch.collect())
	*/

	}

process ENCODE_METADATA {
label "process_single"
tag "${meta.id?:""}"
input:
	val(meta)
output:
	tuple val(meta),path("encode.metadata.tsv"),emit:tsv
	path("versions.yml"),emit:versions
script:
	def url = task.ext.url?:"https://www.encodeproject.org/metadata/?type=Experiment&status=released&assembly=GRCh38&assay_title=total+RNA-seq&files.file_type=bam"
"""
set -o pipefail

curl -L "${url}" |\\
	sed '1s/ /_/g' > encode.metadata.tsv

test -s encode.metadata.tsv

cat << EOF > versions.yml
${task.process}:
    url: "${url}"
	curl: \$( curl --version |head -n1| cut -d ' ' -f2)
EOF
"""
}

process FETCH_GENES {
	label "process_single"
	conda "${moduleDir}/../../conda/bioinfo.01.yml"
	afterScript "rm -rf TMP"
	input:
		tuple val(meta1),path(gtf),path(idx)
		val(meta)
	output:
		tuple val(meta),path("*.gtf.gz"),path("*.gtf.gz.tbi"),emit:gtf
		path("versions.yml"),emit:versions
	script:
		//if(!(meta.contig)) throw new IllegalArgumentException("meta.contig missing");
		//if(!(meta.start)) throw new IllegalArgumentException("meta.start missing");
		//if(!(meta.end)) throw new IllegalArgumentException("meta.end missing");
		def prefix = task.ext.prefix?:"${meta.contig}_${meta.start}_${meta.end}"
	"""
	mkdir -p TMP

	tabix "${gtf}" "${meta.contig}:${meta.start}-${meta.end}" |\\
		awk -F '\t' '(\$3=="gene") {printf("%s\t%s\t%s\\n",\$1,\$4,\$5);}' |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\\
		bedtools merge > TMP/genes.tsv
	
	if ! test -s TMP/genes.tsv
	then
		echo "${meta.contig}\t0\t1" > TMP/genes.bed
	fi

	S=\$(cut -f 2 TMP/genes.tsv | sort -n | head -n 1)
	E=\$(cut -f 3 TMP/genes.tsv | sort -n | tail -n 1)

	tabix "${gtf}" "${meta.contig}:${meta.start}-${meta.end}" |\\
		awk -F '\t' -vS="\${S}" -vE="\${E}" '(\$1=="${meta.config}" && int(\$4)>=int(S) && int(\$5)<=int(E))' |\\
		sort -T TMP -t '\t' -k1,1 -k4,4n |\\
		bgzip > ${prefix}.gtf.gz

	tabix -f -p gff ${prefix}.gtf.gz

cat << EOF > versions.yml
${task.process}:
	tabix: "todo"
EOF
	"""
}

process GET_EXONS {
	label "process_single"
	conda "${moduleDir}/../../conda/bioinfo.01.yml"
	afterScript "rm -rf TMP"
	input:
		tuple val(meta1),path(gtf),path(idx)
		tuple val(meta),path(bed)
	output:
		tuple val(meta),path("*.bed"),emit:bed
		path("versions.yml"),emit:versions
	script:
	"""
	        tabix --regions "${bed}" "${gtf}"  |\\
                awk -F '\t' '(\$3=="exon") {printf("%s\t%s\t%s\\n",\$1,int(\$4)-1,\$5);}' |\\
		sort -T TMP -t '\t' -k1,1 -k4,4n |\\
		bedtools merge > exons.bed

	test -s exons.bed

cat << EOF > versions.yml
${task.process}:
	tabix: "todo"
EOF	"""
	
	}


process SCAN_FOR_SPLICE_EVENTS {
tag "pool ${meta.id}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta),path("BAMS/*")
output:
	tuple val(meta),path("*.bed"),emit:bed
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:"${meta.id}";
"""
mkdir -p TMP

find ./BAMS/ -name "*am" | while read F
do
	jvarkit  -Xmx${task.memory.giga}g bioalcidaejdk --nocode -f "${moduleDir}/extract.code" "\${F}" |\\
		sort -T TMP |\\
		uniq -c |\\
		awk -vS="\${F}" '{printf("%s\t%s\t%s\t%s\t%s\\n",\$2,\$3,\$4,\$1,S);}' |\\
		sed 's%\\./BAMS/%%' | sed 's/\\.bam\$//' >> TMP/jeter.tsv 
done

mv -v TMP/jeter.tsv "${prefix}.bed"

cat << EOF > versions.yml
${task.process}:
	jvarkit: "todo"
EOF
"""
}

process SCAN_FOR_CRYPTIC {
tag "pool ${meta.id}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(tsv)
        tuple val(meta),path("BAMS/*")
output:
        tuple val(meta),path("*.junctions.tsv"),emit:tsv
        path("versions.yml"),emit:versions
script:
        def prefix = task.ext.prefix?:"${meta.id}";
	def undefined="[ATGCNatgcn]{5,20}"
"""
mkdir -p TMP

# bases a droite et a gauche, utilise option -o de grep
awk -F '\t' '{printf("%s_5prime\t${undefined}%s\\n",\$1,\$2);}' "${tsv}" >  TMP/patterns.txt
awk -F '\t' '{printf("%s_3prime\t%s${undefined}\\n",\$1,\$3);}' "${tsv}" >> TMP/patterns.txt
awk -F '\t' '{printf("s|%s|[%s_5prime]|i\\n",\$2,\$1);}' "${tsv}" >  TMP/sed.txt
awk -F '\t' '{printf("s|%s|[%s_3prime]|i\\n",\$3,\$1);}' "${tsv}" >> TMP/sed.txt

set +o pipefail



find ./BAMS/ -name "*am" | while read F
do
	cat TMP/patterns.txt | while read NAME PATTERN
	do
		samtools view -F 3844 "\${F}" |\\
			cut -f 10 |\\
			grep -iEo "\${PATTERN}" |\\
			sed -f TMP/sed.txt |\\
                	sort -T TMP |\\
	                uniq -c |\\
        	        awk -vS="\${F}" -vN="\${NAME}" -vP="\${PATTERN}" '{printf("%s\t%s\t%s\t%s\t%s\\n",\$2,N,P,\$1,S);}' |\\
                	sed 's%\\./BAMS/%%' |\\
			sed 's/\\.bam\$//' >> TMP/jeter.tsv
	done
done

mv -v TMP/jeter.tsv "${prefix}.junctions.tsv"

cat << EOF > versions.yml
${task.process}:
        jvarkit: "todo"
EOF
"""
}





process CAT_SPLICE_EVENTS {
tag "${meta.id?:""}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta),path("BEDS/*")
output:
	tuple val(meta),path("*.bed"),emit:bed
	path("versions.yml"),emit:versions
script:
	def prefix=task.ext.prefix?:"${meta.id}" 
"""
mkdir -p TMP
echo -e "#CHROM\tSTART0\tEND\tCOUNT\tBAM" > TMP/jeter.bed
find ./BEDS/ -name "*.bed" -exec cat '{}' ';' |\\
	sort -t '\t' -T TMP -k1,1 -k2,2n -k3,3n >> TMP/jeter.bed

mv TMP/jeter.bed "${prefix}.bed"

touch versions.yml
"""
}

process CAT_CRYPTIC {tag "${meta.id?:""}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
        tuple val(meta),path("BEDS/*")
output:
        tuple val(meta),path("*.tsv"),emit:bed
        path("versions.yml"),emit:versions
script:
        def prefix=task.ext.prefix?:"${meta.id}"
"""
mkdir -p TMP
echo -e "DNA\tNAME\tDNA\tCOUNT\tBAM" > TMP/jeter.tsv
find ./BEDS/ -name "*.tsv" -exec cat '{}' ';' |\\
        sort -t '\t' -T TMP -k1,1 -k2,2n -k3,3 >> TMP/jeter.tsv

mv TMP/jeter.tsv "${prefix}.tsv"

touch versions.yml
"""
}



process PLOT_SVG {
tag "${meta.id} ${meta.contig}:${meta.start}-${meta.end}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta),path(bam),path(bai),path(gtf),path(gtf_tbi)
output:
	tuple val(meta),path("*.svg.gz", arity: '0..*'),optional:true,emit:svg
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:"${meta.id}";
"""
mkdir -p TMP/OUT

jvarkit  plotsashimi --hyperlink hg38 --skip-empty \\
		-r '${meta.contig}:${meta.start}-${meta.end}' \\
		--gzip \\
		--gtf  "${gtf}" \\
		-o "TMP/OUT" "${bam}"

find TMP/OUT -type f -name "*.svg.gz" | while read F
do
	MD5=`md5sum "\${F}" | awk '{print \$1}'`
	mv -v "\${F}" "./${meta.id}.${meta.contig}_${meta.start}_${meta.end}.\${MD5}.svg.gz"
done

cat << EOF > versions.yml
${task.process}:
	jvarkit: "todo"
EOF
"""
}

