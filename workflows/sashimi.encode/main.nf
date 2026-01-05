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
nextflow.enable.dsl=2



include {runOnComplete               } from '../../modules/utils/functions.nf'
include {JVARKIT_BAM_WITHOUT_BAI     } from '../../modules/jvarkit/bamwithoutbai'
include {BAM_QC                      } from '../../subworkflows/bamqc'
include {MULTIQC                     } from '../../modules/multiqc'
include {COMPILE_VERSIONS            } from '../../modules/versions'
include {BATIK_DOWNLOAD              } from '../../modules/batik/download'
include {BATIK_RASTERIZE             } from '../../modules/batik/rasterize'
include {GHOSTSCRIPT_MERGE           } from '../../modules/gs/merge'
include {PREPARE_REFERENCE           } from '../../subworkflows/samtools/prepare.ref'
runOnComplete(workflow)

workflow {
	version_ch = Channel.empty()
	def hash_ref= [
		id: file(params.fasta).baseName,
		name: file(params.fasta).baseName,
		ucsc_name: (params.ucsc_name?:"undefined"),
		ensembl_name: "GRCh38"
		]
	def fasta = [ hash_ref, file(params.fasta)]
	PREPARE_REFERENCE(hash_ref,Channel.of(fasta))
	fasta = PREPARE_REFERENCE.out.fasta.first()
	version_ch = version_ch.mix(PREPARE_REFERENCE.out.versions)
	
	def gtf = [hash_ref,file(params.gtf),file(params.gtf+".tbi")]
	def bed  = [ hash_ref, file(params.bed)]
	if(!hash_ref.ucsc_name.equals("hg38")) throw new IllegalArgumentException("WANT hg38");
	

	
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
		PREPARE_REFERENCE.out.dict.first(),
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
		PREPARE_REFERENCE.out.fai.first(),
		PREPARE_REFERENCE.out.dict.first(),
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
    MULTIQC(
	[[id:"no_mqc_config"],[]],
	multiqc.map{it[1]}.collect().map{[[id:"sashimi.encode"],it]}
	)

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
stub:
"""
touch versions.yml
cat << EOF | tr "|" "\t" | sed '1s/ /_/g' > encode.metadata.tsv
File accession|File format|File type|File format type|Output type|File assembly|Experiment accession|Assay|Donor(s)|Biosample term id|Biosample term name|Biosample type|Biosample organism|Biosample treatments|Biosample treatments amount|Biosample treatments duration|Biosample genetic modifications methods|Biosample genetic modifications categories|Biosample genetic modifications targets|Biosample genetic modifications gene targets|Biosample genetic modifications site coordinates|Biosample genetic modifications zygosity|Experiment target|Library made from|Library depleted in|Library extraction method|Library lysis method|Library crosslinking method|Library strand specific|Experiment date released|Project|RBNS protein concentration|Library fragmentation method|Library size range|Biological replicate(s)|Technical replicate(s)|Read length|Mapped read length|Run type|Paired end|Paired with|Index of|Derived from|Size|Lab|md5sum|dbxrefs|File download URL|Genome annotation|Platform|Controlled by|File Status|s3_uri|Azure URL|File analysis title|File analysis status|Audit WARNING|Audit NOT_COMPLIANT|Audit ERROR
ENCFF954JBN|bam|bam||alignments|GRCh38|ENCSR000AAA|total RNA-seq|/human-donors/ENCDO018AAA/|CL:0002539|aortic smooth muscle cell|primary cell|Homo sapiens|||||||||||RNA|rRNA||||reverse|2014-06-30|ENCODE||see document|>200|2|2_1||101|||||/files/ENCFF001RAL/, /files/ENCFF001RAM/, /files/ENCFF742NER/|63097172952|ENCODE Processing Pipeline|7b05be865bad92c12c3ea5371aebf6b9||https://www.encodeproject.org/files/ENCFF954JBN/@@download/ENCFF954JBN.bam|V24|||released|s3://encode-public/2016/02/24/8ad79802-f2bf-4a7d-a3c5-db1914a8df83/ENCFF954JBN.bam|https://datasetencode.blob.core.windows.net/dataset/2016/02/24/8ad79802-f2bf-4a7d-a3c5-db1914a8df83/ENCFF954JBN.bam?sv=2019-10-10&si=prod&sr=c&sig=9qSQZo4ggrCNpybBExU8SypuUZV33igI11xw0P7rB3c%3D|ENCODE3 GRCh38 V24|archived|||
ENCFF736LKP|bam|bam||transcriptome alignments|GRCh38|ENCSR000AAA|total RNA-seq|/human-donors/ENCDO018AAA/|CL:0002539|aortic smooth muscle cell|primary cell|Homo sapiens|||||||||||RNA|rRNA||||reverse|2014-06-30|ENCODE||see document|>200|2|2_1||101|||||/files/ENCFF001RAL/, /files/ENCFF001RAM/, /files/ENCFF742NER/|30594465343|ENCODE Processing Pipeline|6c1ffa5580a384e8de8b6332fd7e8fce||https://www.encodeproject.org/files/ENCFF736LKP/@@download/ENCFF736LKP.bam|V24|||released|s3://encode-public/2016/02/24/999e11b8-c8f1-4551-bb68-5759e3a79707/ENCFF736LKP.bam|https://datasetencode.blob.core.windows.net/dataset/2016/02/24/999e11b8-c8f1-4551-bb68-5759e3a79707/ENCFF736LKP.bam?sv=2019-10-10&si=prod&sr=c&sig=9qSQZo4ggrCNpybBExU8SypuUZV33igI11xw0P7rB3c%3D|ENCODE3 GRCh38 V24|archived|||
ENCFF138IAX|bam|bam||alignments|GRCh38|ENCSR000AAA|total RNA-seq|/human-donors/ENCDO018AAA/|CL:0002539|aortic smooth muscle cell|primary cell|Homo sapiens|||||||||||RNA|rRNA||||reverse|2014-06-30|ENCODE||see document|>200|2|2_1||101|||||/files/ENCFF598IDH/, /files/ENCFF001RAL/, /files/ENCFF001RAM/|63099808353|ENCODE Processing Pipeline|072b091ff012a56d3bdc535073863a7f||https://www.encodeproject.org/files/ENCFF138IAX/@@download/ENCFF138IAX.bam|V29|||released|s3://encode-public/2021/04/15/a60cbed2-62bd-4e46-bad2-15b14f46b4e8/ENCFF138IAX.bam|https://datasetencode.blob.core.windows.net/dataset/2021/04/15/a60cbed2-62bd-4e46-bad2-15b14f46b4e8/ENCFF138IAX.bam?sv=2019-10-10&si=prod&sr=c&sig=9qSQZo4ggrCNpybBExU8SypuUZV33igI11xw0P7rB3c%3D|ENCODE4 v1.2.3 GRCh38 V29|released|||

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
stub:
def prefix = task.ext.prefix?:"${meta.contig}_${meta.start}_${meta.end}"
"""
touch versions.yml ${prefix}.gtf.gz ${prefix}.gtf.gz.tbi
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
EOF
"""
stub:
"""
cp ${bed} exons.bed
touch versions.yml
"""	
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
stub:
def prefix = task.ext.prefix?:"${meta.id}"
"""
touch versions.yml "${prefix}.bed"
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

stub:
def prefix = task.ext.prefix?:"${meta.id}"
"""
touch versions.yml "${prefix}.junctions.tsv"
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

stub:
def prefix = task.ext.prefix?:"${meta.id}"
"""
touch versions.yml "${prefix}.bed"
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

stub:
def prefix = task.ext.prefix?:"${meta.id}"
"""
touch versions.yml "${prefix}.tsv"
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

stub:
def prefix = task.ext.prefix?:"${meta.id}"
"""
touch versions.yml "${prefix}.svg.gz"
"""
}

