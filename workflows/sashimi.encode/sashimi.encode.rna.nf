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
nextflow.enable.dsl=2

params.prefix = ""
params.publishDir = ""
params.bed = "NO_FILE"
params.gtf = "NO_FILE"


include {moduleLoad;runOnComplete} from '../../modules/utils/functions.nf'
include {BATIK_DOWNLOAD_01} from '../../modules/batik/batik.download.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {SIMPLE_ZIP_01} from '../../modules/utils/zip.simple.01.nf'
workflow {
	c1_ch = SASHIMI_ENCODE_01(params, Channel.fromPath(params.bed), file(params.gtf))
	html =  VERSION_TO_HTML(params,c1_ch.version)
	zip_ch = Channel.empty().
			mix(c1_ch.pdfs).
			mix(c1_ch.version).
			mix(c1_ch.junctions).
			mix(html.html)
	zip_ch = SIMPLE_ZIP_01([:] ,zip_ch.collect())
	}


runOnComplete(workflow)

workflow SASHIMI_ENCODE_01 {
	take:
		meta
		bed
		gtf
	main:
		version_ch = Channel.empty()
		batik_ch = BATIK_DOWNLOAD_01(meta)
		version_ch = version_ch.mix(batik_ch.version)


		encode_ch = ENCODE_METADATA(meta)
		version_ch = version_ch.mix(encode_ch.version)

		all_bams = encode_ch.output.splitCsv(header:true,sep:'\t').
			filter{T->T.File_format.equals("bam") && T.Output_type.equals("alignments") && T.File_assembly.equals("GRCh38") && !T.File_accession.isEmpty() && !T.File_download_URL.isEmpty()}.
			map{T->[T.File_accession,T.File_download_URL]}

		all_intervals = bed.splitCsv(header: false,sep:'\t',strip:true).
					map{T->[T[0], (T[1] as Integer) +1,T[2] as Integer ]}

		extractor_ch = EXTRACTOR(meta)
		version_ch = version_ch.mix(extractor_ch.version)

		bam_ch = DOWNLOAD_BAM_01(meta, batik_ch.rasterizer, gtf, extractor_ch.output, all_bams.combine(all_intervals))
		version_ch = version_ch.mix(bam_ch.version)

		merge_pdf = MERGE_PDF(meta,bam_ch.pdf.groupTuple())
		version_ch = version_ch.mix(merge_pdf.version)

		junction_ch = COLLECT_ALL_JUNCTIONS(meta,bam_ch.junctions.collect());
		version_ch = version_ch.mix(junction_ch.version)

		version_ch = MERGE_VERSION(meta, "encodeRNA", "sashimi plot encode RNA", version_ch.collect())
	emit:
		version = version_ch
		pdfs = merge_pdf.output
		junctions = junction_ch.output
	}

process ENCODE_METADATA {
executor "local"
input:
	val(meta)
output:
	path("encode.metadata.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def url = meta.url?:"https://www.encodeproject.org/metadata/?type=Experiment&status=released&assembly=GRCh38&assay_title=total+RNA-seq&files.file_type=bam"
"""
set -o pipefail
wget -q -O - "${url}" |\
awk -F '\t' '{OFS="\t";if(NR==1) {for(i=1;i<=NF;i++) gsub(/[ ]/,"_",\$i);} print;}' > encode.metadata.tsv

test -s encode.metadata.tsv

cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">download encode metadata</entry>
        <entry key="url"><url>${url}</url></entry>
</properties>
EOF
"""
}

process EXTRACTOR {
	executor "local"
	input:
		val(meta)
	output:
		path("extract.code"),emit:output
		path("version.xml"),emit:version
	script:
"""

cat << EOF > extract.code
stream().
filter(R->!R.getReadUnmappedFlag()).
forEach(rec->{
	final Cigar cigar=rec.getCigar();
	int refPos1=rec.getAlignmentStart();
	for(int cIdx=0;cIdx< cigar.numCigarElements();++cIdx)
		{
		final CigarElement ce=cigar.getCigarElement(cIdx);
		switch(ce.getOperator())
			{
			case S: break;
			case I: break;
			case N:
				{
				if(cIdx+1<cigar.numCigarElements() &&
					cigar.getCigarElement(cIdx+1).getOperator().isAlignment())
					{
					out.print(rec.getContig());
					out.print("\t");
					out.print(refPos1-1);
					out.print("\t");
					out.print(refPos1+ce.getLength()-1);
					out.println();
					}
				refPos1+=ce.getLength();	
				break;
				}
			case D:
			case X:
			case EQ:
			case M:
				{
				refPos1+=ce.getLength();
				break;
				}
			case H:case P: break;//ignore
			default:throw new IllegalStateException("operator not handled. ops.");
			}
		}
	});

EOF


cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">script extracting junctions locations</entry>
</properties>
EOF

"""
}

process DOWNLOAD_BAM_01 {
        tag "${sample} ${contig}:${start}-${end} ${url}"
        cache 'lenient'
	afterScript 'rm -rf TMP'
	errorStrategy  {task.exitStatus!=0 && task.attempt<3 ? 'retry' : 'ignore' }
	memory "5g"
	maxForks 5
        input:
		val(meta)
		path(rasterizer_jar)
		path(gtf)
		path(code)
                tuple val(sample),val(url),val(contig),val(start),val(end)
        output:
                path("${sample}_${contig}_${start}_${end}.mf"),emit:manifest
                path("${sample}_${contig}_${start}_${end}.tsv"),emit:junctions
		tuple val("${contig}_${start}_${end}"),path("all.pdf.csv"),emit:pdf
                path("version.xml"),emit:version
	when:
		!sample.equals("ENCFF337XBW") //broken bam
        script:
		def proxyHost = "proxy-upgrade.univ-nantes.prive"
		def proxyPort = "3128"
		def proxy = meta.proxy?:" -Dhttp.proxyHost=${proxyHost} -Dhttp.proxyPort=${proxyPort}  -Dhttps.proxyHost=${proxyHost} -Dhttps.proxyPort=${proxyPort} "

        """
	hostname 1>&2
	${moduleLoad("jvarkit samtools tabix")}
	set -o pipefail

	mkdir TMP
	java  -Xmx${task.memory.giga}g ${proxy} -jar \${JVARKIT_DIST}/bamwithoutbai.jar  -r "${contig}:${start}-${end}"  '${url}' |\
	samtools addreplacerg -r '@RG\\tID:${sample}\\tSM:${sample}' -o TMP/jeter.bam -O BAM -
	samtools index TMP/jeter.bam

	tabix "${gtf.toRealPath()}" "${contig}" | java -Xmx${task.memory.giga}g -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "TMP/jeter.bam" --convert SKIP > TMP/tmp.gtf

	# bioalcidae
	java  -Xmx${task.memory.giga}g -jar \${JVARKIT_DIST}/bioalcidaejdk.jar --nocode -f "${code}" TMP/jeter.bam | sort -T TMP |\
		uniq -c |\
		awk '{printf("%s\\t${sample}\\n",\$0);}' > ${sample}_${contig}_${start}_${end}.tsv

	mkdir OUT
	java -jar \${JVARKIT_DIST}/plotsashimi.jar --hyperlink hg38 --skip-empty -r '${contig}:${start}-${end}' \
		--gzip --gtf  TMP/tmp.gtf -m "${sample}_${contig}_${start}_${end}.mf" -o \${PWD}/OUT TMP/jeter.bam

	# apply batik	
	find ./OUT -type f -name "*.svg.gz" | xargs --no-run-if-empty java  -Xmx${task.memory.giga}g -jar ${rasterizer_jar} -m "application/pdf" 

	# find generated pdfs
	find \${PWD}/OUT -type f -name "*.pdf" > all.pdf.csv

	cat <<- EOF > version.xml
	<properties id="${task.process}">
        	<entry key="name">${task.process}</entry>
	        <entry key="description">plot sashimi</entry>
        	<entry key="sample">${sample}</entry>
        	<entry key="url"><url>${url}</url></entry>
        	<entry key="interval">${contig}:${start}-${end}</entry>
        	<entry key="plotsashimi.version">\$(java -jar \${JVARKIT_DIST}/plotsashimi.jar --version)</entry>
	</properties>
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
