params.prefix = ""
params.publishDir = ""
params.bed = "NO_FILE"
params.gtf = "NO_FILE"


workflow {
	SASHIMI_ENCODE_01(params, file(params.bed), file(params.gtf))
	}


workflow SASHIMI_ENCODE_01 {
	take:
		meta
		bed
		gtf
	script:
		version_ch = Channel.empty()

		encode_ch = ENCODE_METADATA(meta)
		version_ch = version_ch.mix(encode_ch.version)

		all_bams = encode_ch.output.splitCsv(header:true,sep:'\t').
			filter{T->T.File_format.equals("bam") && T.Output_type.equals("alignments") && T.File_assembly.equals("GRCh38")}.
			map{T->[T.File_acccession,T.File_download_URL]}

		all_intervals = bed.splitCsv(header: false,sep:'\t',strip:true).
					map{T->[T[0], (T[1] as Integer) +1,T[2] as Integer ]}

		extractor_ch = EXTRACTOR(meta)
		version_ch = version_ch.mix(extractor_ch.version)

		bam_ch = DOWNLOAD_BAM_01(meta, extractor_ch.output, extractor_ch.gtf, all_bams.combine(all_intervals))
		version_ch = version_ch.mix(bam_ch.version)
	emit:
		version = version_ch
	}

process ENCODE_METADATA {
input:
	val(meta)
output:
	path("encode.metadata.tsv"),emit:output
	path("version.xml"),emit:version
script:
"""
set -o pipefail
wget -q -O - "https://www.encodeproject.org/metadata/?type=Experiment&status=released&assembly=GRCh38&assay_title=total+RNA-seq&files.file_type=bam" |\
awk -F '\t' '{OFS="\t";if(NR==1) {for(i=1;i<=NF;i++) gsub(/[ ]/,"_",\$i);} print;}' > encode.metadata.tsv

test -s encode.metadata.tsv

cat << EOF > version.xml
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

"""
}

process DOWNLOAD_BAM_01 {
        tag "${sample} ${contig}:${start}-${end}"
        cache 'lenient'
	afterScript 'rm -rf TMP'
	maxForks 1
        input:
		val(meta)
		path(gtf)
		path(code)
                tuple val(sample),val(url),val(contig),val(start),val(end)
        output:
                path("${sample}_${contig}_${start}_${end}.mf"),emit:manifest
                path("${sample}_${contig}_${start}_${end}.tsv"),emit:manifest
        script:
		def proxyHost = "cache.ha.univ-nantes.fr"
		def proxyPort = "3128"
		def proxy = meta.proxy?:" -Dhttp.proxyHost=${proxyHost} -Dhttp.proxyPort=${proxyPort}  -Dhttps.proxyHost=${proxyHost} -Dhttps.proxyPort=${proxyPort} "

        """
	hostname 1>&2
	${moduleLoad("jvarkit samtools")}
	set -o pipefail

	mkdir TMP
	java ${proxy} -jar \${JVARKIT_DIST}/bamwithoutbai.jar  -r "${contig}:${start}-${end}"  '${url}' |\
	samtools addreplacerg -r '@RG\\tID:${sample}\\tSM:${sample}' -o TMP/jeter.bam -O BAM -
	samtools index TMP/jeter.bam

	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "TMP/jeter.bam" --convert SKIP "${gtf}" > TMP/tmp.gtf

	# bioalcidae
	java -jar \${JVARKIT_DIST}/bioalcidaejdk.jar --nocode -f "${code}" TMP/jeter.bam | sort -T TMP |\
		uniq -c |\
		awk '{printf("%s\\t${sample}\\n",\$0);}' > ${sample}_${contig}_${start}_${end}.tsv

	mkdir OUT
	java -jar \${JVARKIT_DIST}/plotsashimi.jar --hyperlink hg38 --skip-empty -r '${contig}:${start}-${end}' \
		--gzip --gtf  TMP/tmp.gtf -m "${sample}_${contig}_${start}_${end}.mf" -o \${PWD}/OUT TMP/jeter.bam
	
        """
        }

process collectAlljunction {
tag "N=${L.size()}"
input:
	val L from junctions.collect()
output:
	file("${params.prefix}junctions.tsv") into concat_junctions
script:
"""
cat << __EOF__ > jeter.txt
${L.join("\n")}
__EOF__

touch "${params.prefix}junctions.tsv"

cat jeter.txt | while read F
do
	awk '{printf("%s\t%s\t%s\t%s\t%s\\n",\$2,\$3,\$4,\$5,\$1);}' \$F >> "${params.prefix}junctions.tsv"
done

rm jeter.txt
"""
}

process collectPdf {
tag "collect as PDFN=${L.size()}"
publishDir params.publishDir , mode: 'copy', overwrite: true
input:
	val L from manifest.collect()
output:
	file("${params.prefix}all.pdf") into concat_manifest
script:
"""
cat << __EOF__ > jeter.txt
${L.join("\n")}
__EOF__

# le deuxieme cat pour eviter erreur si pas de vcf dans fichier
cat jeter.txt | while read F
do
	grep -v "#" "\$F" | cat >> jeter2.txt
done

# methode recursive parce que ghostview devient lent avec gros pdf
i=0;
sort -t '\t' -k1,1 -k2,2n -k6,6  jeter2.txt | cut -f 7 | while read V
do
	echo \$V 1>&2
	inkscape --without-gui  --file=\$V --export-pdf=chunk.\${i}.pdf
	inkscape --without-gui  --file=\$V --export-width=150 --export-background=white --export-png=chunk.\${i}.png
	((i=i+1))
done



gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=${params.prefix}all.pdf `ls chunk.*.pdf | sort -V` 

rm chunk.*.pdf

"""
}

