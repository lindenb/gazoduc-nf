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
		.filter{!it[1].equals("ENCFF337XBW") }//broken bam
	
	
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
	
	BAM_QC(
		hash_ref,
		fasta,
		fai,
		dict,
		Channel.of(bed).first(),
		all_bams
		)
	version_ch = version_ch.mix(BAM_QC.out.versions)
	multiqc = multiqc.mix(BAM_QC.out.multiqc)


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


process DOWNLOAD_BAM_01 {
    tag "${sample} ${contig}:${start}-${end} ${url}"
        cache 'lenienhostname 1>&2t'
	afterScript 'rm -rf TMP'
	errorStrategy  {task.exitStatus!=0 && task.attempt<3 ? 'retry' : 'ignore' }
	memory "5g"
	maxForks 5
        input:
			path(gtf)
            tuple val(meta),val(sample),val(url),val(contig),val(start),val(end)
        output:
            tuple val(meta),path("*.mf"),emit:manifest
            tuple val(meta),path("*.tsv"),emit:junctions
			tuple val(meta),path("*.svg.gz"),optional:true,emit:svg
			path("versions.yml"),emit:versions
    script:
 """
	hostname 1>&2

	mkdir -p TMP
	jvarkit  -Xmx${task.memory.giga}g \${JAVA_PROXY} bamwithoutbai -r "${contig}:${start}-${end}"  '${url}' |\\
		samtools addreplacerg -r '@RG\\tID:${sample}\\tSM:${sample}' -o TMP/jeter.bam -O BAM -
	
	samtools index TMP/jeter.bam

	tabix "${gtf.toRealPath()}" "${contig}" |\\
		jvarkit -Xmx${task.memory.giga}g bedrenamechr -f "TMP/jeter.bam" --convert SKIP > TMP/tmp.gtf

	# bioalcidae
	jvarkit  -Xmx${task.memory.giga}g bioalcidaejdk --nocode -f "${moduleDir}/extract.code" TMP/jeter.bam |\\
		sort -T TMP |\\
		uniq -c |\\
		awk '{printf("%s\\t${sample}\\n",\$0);}' > ${sample}_${contig}_${start}_${end}.tsv

	mkdir -p OUT
	jvarkit  plotsashimi --hyperlink hg38 --skip-empty -r '${contig}:${start}-${end}' \\
		--gzip --gtf  TMP/tmp.gtf -m "${sample}_${contig}_${start}_${end}.mf" -o \${PWD}/OUT TMP/jeter.bam

	# apply batik	
	#find ./OUT -type f -name "*.svg.gz" | xargs --no-run-if-empty java  -Xmx${task.memory.giga}g -jar ${rasterizer_jar} -m "application/pdf" 

	# find generated pdfs
	#find \${PWD}/OUT -type f -name "*.pdf" > all.pdf.csv

cat <<- EOF > versions.yml
"${task.process}"
	jvarkit: todo
EOF
"""
}

process COLLECT_ALL_JUNCTIONS {
tag "N=${L.size()}"
input:
	val(meta)
	val(L)
output:
	path("${params.prefix}junctions.tsv.gz"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
set -o pipefail


cat << __EOF__ > jeter.txt
${L.join("\n")}
__EOF__

touch "${params.prefix}junctions.tsv"

cat jeter.txt | while read F
do
	awk '{printf("%s\t%s\t%s\t%s\t%s\\n",\$2,\$3,\$4,\$5,\$1);}' \$F >> "${params.prefix}junctions.tsv"
done

gzip --best "${params.prefix}junctions.tsv"

rm jeter.txt
	
cat <<- EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
	<entry key="description">merge junctions</entry>
        <entry key="count">${L.size()}</entry>
</properties>
EOF
"""
}

process MERGE_PDF {
tag "${interval} N=${L.size()}"
input:
	val(meta)
	tuple val(interval),val(L)
output:
	path("${params.prefix?:""}${interval}.pdf"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
set -o pipefail

cat ${L.join(" ")} | awk -F '/' '{printf("%s,%s\\n",\$NF,\$0);}' |\
	sort -t, -k1,1 -T. | cut -d, -f2  > jeter.list
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile="${params.prefix?:""}${interval}.pdf" @jeter.list
rm jeter.list

cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">merge pdfs</entry>
        <entry key="gs.version">\$(gs --version)</entry>
</properties>
EOF
"""
}
