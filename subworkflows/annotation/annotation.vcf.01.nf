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

include {getVersionCmd;isBlank;moduleLoad} from '../../modules/utils/functions.nf'
include {CONCAT_FILES_01} from '../../modules/utils/concat.files.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_ANNOTATE_SOURCES} from './bcftools.annotate.sources.01.nf'

boolean hasFeature(key) {
	def with = "with_"+key
	if(!params.containsKey("annotations")) return false;
	if(!params.annotations.containsKey(with)) {
		log.warn("undefined params.annotations."+with);
		}
	return (params.annotations[with] as boolean)
	}

workflow ANNOTATE_VCF_01 {
	take:
		meta
		genomeId
		rows
	main:
		version_ch = Channel.empty()


		annotate_src_ch = BCFTOOLS_ANNOTATE_SOURCES([:], genomeId, file("NO_FILE"))
		
		rows = rows.combine(annotate_src_ch.annotations_ch.toList().map{T->[annotate:T]}).
			map{T->T[0].plus(T[1])}
		
	
		version_ch = version_ch.mix(annotate_src_ch.version)

		if(params.genomes[genomeId].containsKey("gtf")) {

			if(hasFeature("snpeff") && !hasFeature("native_snpeff") && params.genomes[genomeId].containsKey("gtf")) {
				snpeff_db = BUILD_SNPEFF([:],genomeId)
        		        version_ch = version_ch.mix(snpeff_db.version)

				rows_ch = rows.combine(snpeff_db.output).
					map{T->T[0].plus("snpeff_data":T[1])}
				}
			}

		annot_vcf_ch = APPLY_ANNOTATION(
			[:],
			genomeId,
			rows
			)

/*



		version_ch = version_ch.mix(annot_vcf_ch.version.first())

		to_file_ch = COLLECT_TO_FILE_01([:],annot_vcf_ch.bedvcf.map{T->T[0]+"\t"+T[1].toRealPath()}.collect())
		version_ch = version_ch.mix(to_file_ch.version)

*/
		version_ch = MERGE_VERSION(meta, "annotation", "VCF annotation", version_ch.collect())
	emit:
		version= version_ch
		//bedvcf = to_file_ch.output
	}


process complementCaptures {
memory "3g"
input:
	val(meta)
output:
	path("norm.captures.tsv"),emit:output
script:
if(isBlank(meta.captures))
"""
touch norm.captures.tsv
"""
else
"""
hostname 1>&2
${moduleLoad("bedtools jvarkit htslib")}
set -o pipefail

touch norm.captures.tsv

cut -f1,2 "${meta.reference}.fai" |\
	sort -t '\t' -k1,1 -k2,2n -T .  > jeter.genome

grep -v "^#" "${meta.captures}" | while read CN CF
	do

		(gunzip -c "\${CF}" || cat "\${CF}") | grep -v -E '^(browser|track|#)' | cut -f1,2,3 |\
		java -Xmx${task.memory.giga}g -Djava.io.tmpdir=.  -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${meta.reference}" --column 1 --convert SKIP  |\
               		sort -t '\t' -k1,1 -k2,2n -T . |\
			bedtools merge |\
			bedtools complement -i - -g jeter.genome |\
			awk -F '\t' '(int(\$2) < int(\$3)) {printf("%s\t1\\n",\$0);}' |\
			bedtools sort -faidx "${meta.reference}.fai" |\
               		bgzip > "\${CN}.bed.gz"

			tabix -f -p bed "\${CN}.bed.gz"

			echo "\${CN}\t\${PWD}/\${CN}.bed.gz" >> norm.captures.tsv
	done

rm jeter.genome
"""
}


/** build snpeff Database from gtf */
process BUILD_SNPEFF {
afterScript "rm -f org.tsv genes.tsv"
memory "10g"
input:
        val(meta)
        val(genomeId)
output:
	path("data"),emit:output
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def dbName = genomeId
	def reference = genome.fasta
	def gtf = genome.gtf
"""
hostname 1>&2
${moduleLoad("snpEff")}
set -o pipefail

mkdir -p "data/${dbName}"
ln -s "${reference}" "data/${dbName}/sequences.fa"

cp "${gtf}"  data/${dbName}/genes.gtf.gz
gunzip data/${dbName}/genes.gtf.gz

# write snpEff contig
cat << EOF > snpEff.config
data.dir = \${PWD}/data/
${dbName}.genome = Human
${dbName}.reference = ${reference}
EOF

# build database
java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar \${SNPEFF_JAR}  build -gtf22 -v "${dbName}" 2> snpeff.errors

rm snpeff.errors

test -s "data/${dbName}/snpEffectPredictor.bin"

rm  data/${dbName}/genes.gtf

cat <<- EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">build custom snpeff from gtf</entry>
        <entry key="gtf">${gtf}</entry>
        <entry key="fasta">${reference}</entry>
        <entry key="output">\${PWD}/snpEff.config</entry>
        <entry key="snpeff">\$(java -jar \${SNPEFF_JAR} -version)</entry>
</properties>
EOF
"""
}


process DOWNLOAD_TOMMO {
afterScript "rm -f tmp.bed"
input:
	val(meta)
	val(reference)
      	path bed from merged_bed
output:
       	path("tommo.bed.gz") into tommo_vcf
        path("tommo.bed.gz.tbi")
script:
if(isHg19(reference) && hasFeature("tommo"))
"""
hostname 1>&2
module load htslib/0.0.0 jvarkit bedtools/0.0.0 bcftools/0.0.0

for U in "tommo-14KJPN-20211208-GRCh37_lifted_from_GRCh38-af-autosome" 
do
	wget -O - "https://jmorp.megabank.tohoku.ac.jp/dj1/datasets/tommo-14kjpn-20211208-af_snvindelall/files/\${U}.vcf.gz?download=true" |\
		gunzip -c|\
		java -jar \${JVARKIT_DIST}/vcfsetdict.jar -R "${meta.reference}" |\
		bcftools norm -m '-' --targets-file "${bed}" -O u |\
		bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\n' >> tmp.bed
done
LC_ALL=C sort -t '\t' -k1,1 -k2,2n -T . tmp.bed | uniq |  bgzip > tommo.bed.gz
tabix -s 1 -b 2 -e 2 tommo.bed.gz
"""
else
"""
module load htslib/0.0.0
touch tommo.bed
bgzip tommo.bed
tabix -s 1 -b 2 -e 2 tommo.bed.gz
"""
}





process APPLY_ANNOTATION {
tag "${row.vcf} ${row.interval}"
afterScript "rm -rf TMP"
memory '5 g'
maxForks 30

input:
	val(meta)
	val(genomeId)
	val(row)
output:
	tuple val(row),path("contig.bcf"),emit:bedvcf
	path("contig.bcf.csi"),emit:index
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
	def vcf=row.vcf
	def pedigree = row.pedigree?:file("NO_FILE")
	def extraBcfTools=""
"""
	hostname 1>&2
	${moduleLoad("bcftools jvarkit snpEff bedtools htslib")}
	set -o pipefail
	set -x
	mkdir -p TMP
	
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Annotation of a VCF file</entry>
		<entry key="bed">${row.bed?:""}</entry>
		<entry key="interval">${row.interval?:""}</entry>
		<entry key="vcf">${vcf}</entry>
		<entry key="bcftools.version">\$(bcftools --version | head -n 2 | paste -sd ' ')</entry>
		<entry key="steps">
	EOF

	# save bed or interval
	if ${row.containsKey("bed")} ; then
		ln -s "${row.bed}" TMP/tmp.bed
	else
		echo "${row.interval}" | awk -F '[:-]' '{printf("%s\t%d\t%s\\n",\$1,int(\$2)-1,\$3);}' > TMP/tmp.bed
	fi
	

	# normalize bed, bcftools doesn't like more than 3 columns...
	cut -f1,2,3 "TMP/tmp.bed" | bedtools sort -faidx "${reference}.fai" > TMP/jeter.123.bed

	# extract variants
	bcftools view ${params.annotations.extraBcfTools} --regions-file "TMP/jeter.123.bed" -O v "${vcf}" -o TMP/jeter1.vcf

	# samples in pedigree
	if ${!pedigree.name.equals("NO_FILE") && hasFeature("keepSamplesInPed")} ; then
		cut -f 2 "${pedigree}"  | sort -T TMP | uniq > TMP/samples.a
		bcftools query -l TMP/jeter1.vcf  |	sort -T TMP | uniq > TMP/samples.b
		comm -12 TMP/samples.a TMP/samples.b > TMP/samples.c
		test -s TMP/samples.c
		bcftools view --samples-file TMP/samples.c -O u TMP/jeter1.vcf |\
			bcftools view -i 'AC[*]>0' > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi

	# annotate variant ID
	if ${hasFeature("dbsnp_id") && genome.containsKey("dbsnp")} ; then

		# file must be indexed
		bcftools view -O b -o TMP/jeter1.bcf TMP/jeter1.vcf
		bcftools index TMP/jeter1.bcf
		
		bcftools annotate --annotations "${genome.dbsnp}"  --regions-file "TMP/jeter.123.bed"  -c ID -O v -o TMP/jeter1.vcf TMP/jeter1.bcf
		rm TMP/jeter1.bcf TMP/jeter1.bcf.csi
	fi

	# bcftools CSQ
        if ${hasFeature("bcftools_csq") && genome.gff3}  ; then

  		bcftools csq -O v --force --local-csq --ncsq 10000 --fasta-ref "${reference}" --gff-annot "${genome.gff3}" TMP/jeter1.vcf > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi

	# annotations with bcftools annotate
	if ${row.annotate && !row.annotate.isEmpty()} ; then
		#
		${row.annotate.collect{T->"#####\nbcftools annotate --annotations '${T.tabix}' --columns '${T.columns}' --header-lines '${T.header}' --merge-logic '${T.merge_logic}' ${T.minoverlap.equals(".")?"":"--min-overlap '${T.minoverlap}'"} -O v -o TMP/jeter2.vcf TMP/jeter1.vcf\nmv -v TMP/jeter2.vcf TMP/jeter1.vcf"}.join("\n")}

	fi


bcftools sort --max-mem "${task.memory.giga}g" -T TMP  -O b -o contig.bcf TMP/jeter1.vcf
bcftools index contig.bcf


cat <<- EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">annot vcf</entry>
</properties>
EOF
"""
}

// 		


process TODOO {

script:
"""

	# select with jvarkit
	if [ ! -z "${isBlank(meta.select)?"":"Y"}" ] ; then
		java  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/vcffilterjdk.jar --nocode  ${meta.select} TMP/jeter1.vcf > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi
	
	


	# polyx
	if [ ! -z "${hasFeature("polyx") && (meta.polyx as Integer)>0 ?"Y":""}" ] ; then
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar ${jvarkit("vcfpolyx")} --filter ${meta.polyx} --reference ${reference} --tag POLYX TMP/jeter1.vcf > TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi
	

	# allelic ratio
	if [ ! -z "${hasFeature("AD")?"Y":""}" ] ; then
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar ${jvarkit("vcffilterjdk")} \
		--nocode \
		${isSoftFilter(meta,"AD_RATIO")?"--filter HET_BAD_AD_RATIO ":""} \
		-e 'return variant.getGenotypes().stream().filter(G->G.isHet() && G.hasAD() && G.getAD().length==2).allMatch(G->{int array[]=G.getAD();double r= array[1]/(double)(array[0]+array[1]);return (r>=0.2 && r<=0.8) ;});' TMP/jeter1.vcf > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi

	# GQ
	if [ ! -z "${hasFeature("GQ")?"Y":""}" ] ; then
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar ${jvarkit("vcffilterjdk")} \
		--nocode \
		${isSoftFilter(meta,"LOW_GQ")?"--filter LOW_GQ${meta.lowGQ} ":""} \
		-e 'return variant.getGenotypes().stream().filter(G->(G.isHet() || G.isHomVar()) && G.hasGQ()).allMatch(G->G.getGQ()>=${meta.lowGQ});' TMP/jeter1.vcf > TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi

	# DP
	if [ ! -z "${hasFeature("DP")?"Y":""}" ] ; then
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar ${jvarkit("vcffilterjdk")} \
		--nocode \
		${isSoftFilter(meta,"LOW_DP")?"--filter LOW_DP${meta.lowDP} ":""} \
		-e 'return variant.getGenotypes().stream().filter(G->(G.isHet() || G.isHomVar()) && G.hasDP()).allMatch(G->G.getDP()>= ${meta.lowDP} );' TMP/jeter1.vcf > TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi

	# MQ
	if [ ! -z "${hasFeature("MQ")?"Y":""}" ] ; then
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar ${jvarkit("vcffilterjdk")} \
		--nocode \
		${isSoftFilter(meta,"LOW_MQ")?"--filter LOW_MQ${meta.lowMQ} ":""} \
		-e 'return !variant.hasAttribute("MQ") || variant.getAttributeAsDouble("MQ",1000.0) >= ${meta.lowMQ};' TMP/jeter1.vcf > TMP/jeter2.vcf

	mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi

	# MQRankSum
	if [ ! -z "${hasFeature("MQRankSum")?"Y":""}" ] ; then
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar ${jvarkit("vcffilterjdk")} \
		--nocode \
		${isSoftFilter(meta,"BAD_MQRankSum")?"--filter BAD_MQRankSum${meta.mqRankSum} ":""} \
		-e 'return !variant.hasAttribute("MQRankSum") || Math.abs(variant.getAttributeAsDouble("MQRankSum",0.0)) <= ${meta.mqRankSum};' TMP/jeter1.vcf > TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi

	# ReadPosRankSum
	if [ ! -z "${hasFeature("ReadPosRankSum")?"Y":""}" ] ; then
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar ${jvarkit("vcffilterjdk")} \
		--nocode \
		${isSoftFilter(meta,"BAD_ReadPosRankSum")?"--filter BAD_ReadPosRankSum${meta.readPosRankSum} ":""} \
		-e 'return !variant.hasAttribute("ReadPosRankSum") || Math.abs(variant.getAttributeAsDouble("ReadPosRankSum",0.0)) <= ${meta.readPosRankSum};' TMP/jeter1.vcf > TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi

	# SOR
	if [ ! -z "${hasFeature("SOR")?"Y":""}" ] ; then
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar ${jvarkit("vcffilterjdk")} \
		--nocode \
		${isSoftFilter(meta,"BAD_SOR")?"--filter BAD_SOR${meta.sor} ":""} \
		-e 'return !variant.hasAttribute("SOR") || variant.getAttributeAsDouble("SOR",10000.0) <= ${meta.sor};' TMP/jeter1.vcf > TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi


	# ANNOTATIONS
	cat "${annotations_files}" | while read DATABASE COLS HEADER
	do
		echo "#ANNOT \${DATABASE}" 1>&2
		test -f "\${DATABASE}"
		test -f "\${HEADER}"
		test -f "\${DATABASE}.tbi"

		bcftools annotate --annotations "\${DATABASE}" \
			-h "\${HEADER}" \
			-c "\${COLS}" \
			`echo "\${COLS}" | tr "," "\\n" | awk '(\$1!="CHROM" && \$1!="FROM" && \$1!="POS" && \$1="TO" && \$1!="END" && \$1!="." && \$1!="-") {S="unique"; printf(" --merge-logic %s:%s ",\$1,S);}'  ` \
			-O v -o TMP/jeter2.vcf TMP/jeter1.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf

		cat <<- EOF >> version.xml
		<properties>
			<entry key="description">annotate</entry>
			<entry key="header">\${HEADER}</entry>
			<entry key="info.columns">\${COLS}</entry>
		</properties>
		EOF
	done




	if [ ! -z "${isHg19(reference) && hasFeature("snpeff")?"Y":""}"  ] ; then
		java  -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar "\${SNPEFF_JAR}" eff \
			-config  "\${SNPEFF_CONFIG}" \
			-nodownload -noNextProt -noMotif -noInteraction -noLog -noStats -chr chr -i vcf -o vcf "${isHg19(reference)?"GRCh37.75":"TODO"}" TMP/jeter1.vcf > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf


		cat <<- EOF >> version.xml
		<properties>
			<entry key="description">annotation with SnpEFF</entry>
		</properties>
		EOF


	fi


        if [ ! -z "${isHg19(reference) && hasFeature("vep")?"Y":""}" ] ; then
		module load ensembl-vep/104.3
		vep --cache  --format vcf --force_overwrite --output_file STDOUT --no_stats --offline  --dir_cache /LAB-DATA/BiRD/resources/apps/vep  --species homo_sapiens --cache_version 104 --assembly GRCh37 --merged --fasta "${reference}" --use_given_ref --vcf < TMP/jeter1.vcf > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf
		module unload ensembl-vep/104.3


		cat <<- EOF >> version.xml
		<properties>
			<entry key="description">annotation with VEP</entry>
			<entry key="vep.version">104.3</entry>
		</properties>
		EOF


	fi

        if [ ! -z "${(hasFeature("snpeff") || hasFeature("vep")) && !isBlank(meta.soacn) ?"Y":""}" ] ; then
		java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar ${jvarkit("vcffilterso")} \
			${isSoftFilter(meta,"BAD_SO")?"--filterout  BAD_SO":""} \
			--acn "${meta.soacn}" TMP/jeter1.vcf > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf

		cat <<- EOF >> version.xml
		<properties>
			<entry key="description">filtration consequences</entry>
			<entry key="so.acn">${meta.soacn}</entry>
			<entry key="soft.filter">${isSoftFilter(meta,"BAD_SO")}</entry>
			<entry key="vcffilterso.version">\$(java -jar ${jvarkit("vcffilterso")} --version)</entry>
		</properties>
		EOF


	fi

	## cadd
	if [ ! -z  "${isHg19(reference) && hasFeature("cadd")?"Y":""}" ] ; then
        	java  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar ${jvarkit("vcfcadd")} \
			--tabix "/LAB-DATA/BiRD/resources/species/human/krishna.gs.washington.edu/download/CADD/v1.6/whole_genome_SNVs.tsv.gz" TMP/jeter1.vcf > TMP/jeter2.vcf
	      	mv TMP/jeter2.vcf TMP/jeter1.vcf


		cat <<- EOF >> version.xml
		<properties>
			<entry key="description">annotation CADD</entry>
			<entry key="cadd.file">/LAB-DATA/BiRD/resources/species/human/krishna.gs.washington.edu/download/CADD/v1.6/whole_genome_SNVs.tsv.gz</entry>
			<entry key="vcfcadd.version">\$(java -jar ${jvarkit("vcfcadd")} --version)</entry>
		</properties>
		EOF


	fi

	if [ ! -z "${hasFeature("gnomadGenome")?"Y":""}" ] ; then
		java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar ${jvarkit("vcfgnomad")} \
			--bufferSize 10000 \
			--max-af ${meta.gnomadAF} \
			--prefix "${isSoftFilter(meta,"GNOMAD")?"GNOMAD":""}" \
			--gnomad "${getGnomadGenomePath(meta,reference)}" --fields "${meta.gnomadPop}" TMP/jeter1.vcf > TMP/jeter2.vcf
                mv TMP/jeter2.vcf TMP/jeter1.vcf

		cat <<- EOF >> version.xml
		<properties>
			<entry key="description">annotation gnomad genome</entry>
			<entry key="max-AF">${meta.gnomadAF}</entry>
			<entry key="population">${meta.gnomadPop}</entry>
			<entry key="soft.filter">${isSoftFilter(meta,"GNOMAD")}</entry>
			<entry key="gnomad.path">${getGnomadGenomePath(meta,reference)}</entry>
			<entry key="vcfgnomad.version">\$(java -jar ${jvarkit("vcfgnomad")} --version)</entry>
		</properties>
		EOF

	fi

	if [ ! -z "${hasFeature("gnomadExome")?"Y":""}" ] ; then
		java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar ${jvarkit("vcfgnomad")} \
			--bufferSize 10000 \
			--max-af ${meta.gnomadAF} \
			--prefix "${isSoftFilter(meta,"GNOMAD")?"GNOMAD":""}" \
			--gnomad "${getGnomadExomePath(meta,reference)}" --fields "${meta.gnomadPop}" TMP/jeter1.vcf > TMP/jeter2.vcf
                mv TMP/jeter2.vcf TMP/jeter1.vcf


		cat <<- EOF >> version.xml
		<properties>
			<entry key="description">annotation gnomad exome</entry>
			<entry key="max-AF">${meta.gnomadAF}</entry>
			<entry key="population">${meta.gnomadPop}</entry>
			<entry key="soft.filter">${isSoftFilter(meta,"GNOMAD")}</entry>
			<entry key="gnomad.path">${getGnomadExomePath(meta,reference)}</entry>
			<entry key="vcfgnomad.version">\$(java -jar ${jvarkit("vcfgnomad")} --version)</entry>
		</properties>
		EOF


	fi


	# de novo
	if [ ! -z "${pedigree.name.equals("NO_FILE")?"":"Y"}" ] ; then
		# check parents
		awk -F '\t' '(\$3!="0" || \$4!="0")' "${pedigree}" > TMP/parents.txt
		
		if [ -s TMP/parents.txt ] ; then
			java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar -jar \${JVARKIT_DIST}/vcftrio.jar \
				--pedigree "${pedigree}" TMP/jeter1.vcf > TMP/jeter2.vcf
			mv TMP/jeter2.vcf TMP/jeter1.vcf

			cat <<- EOF >> version.xml
			<properties>
				<entry key="description">annotation of trios</entry>
				<entry key="pedigree">${pedigree}</entry>
				<entry key="vcftrio.version">\$(java -jar \${JVARKIT_DIST}/vcftrio.jar --version)</entry>
			</properties>
			EOF

		fi
		rm TMP/parents.txt

		if [ ! -z "${hasFeature("contrast")?"Y":""}" ] ; then
			# check cases and controls
			awk -F '\t' '(\$6=="case" || \$6=="affected" || \$6=="1")' "${pedigree}" | cut -f 2 > TMP/cases.list
			awk -F '\t' '(\$6=="control" || \$6=="unaffected" || \$6=="0")' "${pedigree}" | cut -f 2 > TMP/ctrl.list

			if [ ! -s TMP/cases.list ] && [ ! -s TMP/ctrl.list ] ; then

    				bcftools +contrast -0 TMP/ctrl.list -1 TMP/cases.list \
      					-a PASSOC,FASSOC,NASSOC,NOVELAL,NOVELGT -O b -o TMP/jeter2.vcf TMP/jeter1.vcf
				mv TMP/jeter2.vcf TMP/jeter1.vcf


				cat <<- EOF >> version.xml
				<properties>
					<entry key="description">bcftools contrast</entry>
					<entry key="pedigree">${pedigree}</entry>
					<entry key="bcftools.version">\$(bcftools version | head -n2 | paste -d ' ' -s)</entry>
				</properties>
				EOF

			fi
			rm cases.list ctrl.list
		fi
	fi


	# add external info AF,AC,AN e.g. amalgamion
	if ${hasFeature("joinvcf") && !getAmalgamionVcf(meta,reference).isEmpty()} ; then
		# get interval for this vcf
		bcftools query -f '%CHROM\t%POS0\t%END\\n' TMP/jeter1.vcf |\
			sort -T TMP -t '\t' -k1,1 -k2,2n |\
			bedtools merge > TMP/jeter.bed

		# if bed is empty
		if ! [ -s  TMP/jeter.bed ] ; then
			head -n1 "TMP/jeter.123.bed" | awk -F '\t' '{printf("%s\t0\t1\\n",\$1);}'  > TMP/jeter.bed
		fi

		bcftools norm --regions-file TMP/jeter.bed -O u -m- "${getAmalgamionVcf(meta,reference)}" |\
			bcftools  +fill-tags -O u  -- -t AN,AC,AF |\
			bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\\n' > TMP/jeter.tab
		rm TMP/jeter.bed
		bgzip TMP/jeter.tab
		tabix -s1 -b2 -e2 TMP/jeter.tab.gz

		echo -e '##INFO=<ID=AC${getAmalgamionSuffix(meta)},Number=A,Type=Integer,Description="AC from ${getAmalgamionVcf(meta,reference)}">' > TMP/hdr.txt
		echo -e '##INFO=<ID=AF${getAmalgamionSuffix(meta)},Number=A,Type=Float,Description="AF from ${getAmalgamionVcf(meta,reference)}">' >> TMP/hdr.txt
		echo -e '##INFO=<ID=AN${getAmalgamionSuffix(meta)},Number=1,Type=Integer,Description="AF from ${getAmalgamionVcf(meta,reference)}">' >> TMP/hdr.txt
		
		bcftools annotate  --mark-sites +IN${getAmalgamionSuffix(meta)} -a TMP/jeter.tab.gz -h TMP/hdr.txt -c 'CHROM,POS,REF,ALT,AC${getAmalgamionSuffix(meta)},AN${getAmalgamionSuffix(meta)},AF${getAmalgamionSuffix(meta)}' TMP/jeter1.vcf > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf

		rm TMP/hdr.txt TMP/jeter.tab.gz TMP/jeter.tab.gz.tbi


		cat <<- EOF >> version.xml
		<properties>
			<entry key="description">annotation AC,AF,AN with external VCF</entry>
			<entry key="vcf">${getAmalgamionVcf(meta,reference)}</entry>
			<entry key="bcftools.version">\$(bcftools version | head -n2 | paste -d ' ' -s)</entry>
		</properties>
		EOF


	fi



	bcftools sort --max-mem "${task.memory.giga}g" -T TMP -O b -o TMP/contig.bcf TMP/jeter1.vcf
	bcftools index TMP/contig.bcf

	mv TMP/contig.bcf ./
	mv TMP/contig.bcf.csi ./

	cat <<- EOF >> version.xml
	</entry>
	</properties>
	EOF


"""
stub:
"""
	# tommo VCF japan variants
	if [ "${isHg19(refrence) && hasFeature("tommo")?"Y":"N"}" == "Y" ] ; then

		echo -e '##INFO=<ID=AC_TOMMO,Number=A,Type=Integer,Description="AC from Tommo DB">' > hdr.txt
		echo -e '##INFO=<ID=AF_TOMMO,Number=A,Type=Float,Description="AF from Tommo DB">' >> hdr.txt
		echo -e '##INFO=<ID=AN_TOMMO,Number=1,Type=Integer,Description="AF from Tommo DB">' >> hdr.txt
		
		bcftools annotate  --mark-sites +IN_TOMMO -a "${tommo}" -h hdr.txt -c 'CHROM,POS,REF,ALT,AC_TOMMO,AN_TOMMO,AF_TOMMO' jeter1.vcf > jeter2.vcf
		mv jeter2.vcf jeter1.vcf

		rm hdr.txt
	fi



	# merge external vcf e.g. amalgamion
	if [ ! -z "${isBlank(mergevcf)?"":"Y"}" ] && [ ! -z "${hasFeature("amalgamion")?"Y":""}" ] ; then
		# get interval for this vcf
		bcftools query -f '%CHROM\t%POS0\t%END\\n' jeter1.vcf |\
			sort -T . -t '\t' -k1,1 -k2,2n |\
			bedtools merge > jeter.bed

		# if bed is empty
		if ! [ -s  jeter.bed ] ; then
			head -n1 "jeter.123.bed" | awk -F '\t' '{printf("%s\t0\t1\\n",\$1);}'  > jeter.bed
		fi

		# exclude common samples
		bcftools query -l jeter1.vcf |sort | uniq > samples1.txt
		bcftools query -l "${mergevcf}" |sort | uniq > samples2.txt
		comm -13 samples1.txt samples2.txt > samples3.txt
		test -s samples3.txt
		
		if [ -s  samples1.txt ] ; then
			awk '{printf(" %s variant.getGenotype(\\"%s\\").isNoCall()", (NR==1?" return !(" :" && "),\$0);} END {printf(") ;\\n");} ' samples1.txt > jeter.code
		else
			echo "return true;" > jeter.code
		fi

		rm samples1.txt samples2.txt 

		# get variants to join in that region
		bcftools annotate -x 'ID,QUAL,FILTER,INFO' --regions-file jeter.bed -O u "${mergevcf}"  |\
			bcftools view --samples-file samples3.txt -O b -o jeter4.bcf
		rm samples3.txt jeter.bed
		bcftools index jeter4.bcf
	
		# bcftools merge wants indexed file
		bcftools view -O b -o jeter1.vcf.gz jeter1.vcf
		bcftools index jeter1.vcf.gz

		# run merge
		bcftools merge  --regions-file "jeter.123.bed" -O v -o jeter1.vcf jeter1.vcf.gz jeter4.bcf

		# keep where a genotype of original VCF has one called genotype
		java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar ${jvarkit("vcffilterjdk")} -f jeter.code jeter1.vcf > jeter2.vcf
		mv jeter2.vcf jeter1.vcf
		

		bcftools query -f '%CHROM:%POS\\n' jeter1.vcf  | sort | uniq > tmp1.txt
		bcftools query -f '%CHROM:%POS\\n' jeter1.vcf.gz  | sort | uniq > tmp2.txt
		comm -3 tmp1.txt tmp2.txt >&2
		rm tmp1.txt tmp2.txt
		rm jeter4.bcf jeter4.bcf.csi jeter1.vcf.gz jeter1.vcf.gz.csi
		

	fi



	
	bcftools sort --max-mem "${task.memory.giga}g" -T . -O b -o contig.bcf jeter1.vcf
	bcftools index contig.bcf
	rm jeter1.vcf
"""
}

