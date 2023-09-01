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

include {getVersionCmd;moduleLoad} from '../../modules/utils/functions.nf'
include {CONCAT_FILES_01} from '../../modules/utils/concat.files.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_ANNOTATE_SOURCES} from './bcftools.annotate.sources.01.nf'
include {UTR_ANNOTATOR_DOWNLOAD_01} from '../../modules/vep/utrannotator.download.01.nf'
include {WGET_01 as WGET_PHYLOP} from '../../modules/utils/wget.01.nf'
include {WGET_01 as WGET_PHASTCONS} from '../../modules/utils/wget.01.nf'

boolean isBlank(hash,key) {
	if(hash==null) throw new IllegalArgumentException("[isBlank] hash is null");
	if(!hash.containsKey(key)) return true;
	def v  = hash.get(key);
	if(v==null || v.toString().trim().isEmpty()) return true;
	return false;
	}

boolean hasFeature(key) {
	def with = "with_"+key
	if(!params.containsKey("annotations")) return false;
	if(!params.annotations.containsKey(with)) {
		log.warn("params.annotations."+with+" missing. return false;");
		return false;
		}
	if(isBlank(params.annotations,with)) return false;
	return (params.annotations[with] as boolean)
	}

boolean isSoftFilter(key) {
	def f = "softfilter_"+key
	if(!params.containsKey("annotations")) return true;
	if(!params.annotations.containsKey(f)) {
		log.warn("params.annotations."+f+" missing. return true;");
		return true;
		}
	if(isBlank(params.annotations,f)) return true;
	return (params.annotations[f] as boolean)
	}

workflow ANNOTATE_VCF_01 {
	take:
		meta
		genomeId
		rows
	main:
		version_ch = Channel.empty()


		annotate_src_ch = BCFTOOLS_ANNOTATE_SOURCES([:], genomeId, file("NO_FILE"))
		version_ch = version_ch.mix(annotate_src_ch.version)


		rows = rows.combine(annotate_src_ch.annotations_ch.toSortedList(/* PREDICTIVE ORDER = PREVENT CACHE INVALIDATION */{a,b -> a.toString().compareTo(b.toString())}).map{T->[annotate:T]}).
			map{T->T[0].plus(T[1])}
		
	
		if(hasFeature("snpeff") && !hasFeature("native_snpeff") && !isBlank(params.genomes[genomeId],"gtf")) {
			snpeff_db = BUILD_SNPEFF([:],genomeId)
        		version_ch = version_ch.mix(snpeff_db.version)

			rows = rows.combine(snpeff_db.output).
					map{T->T[0].plus("snpeff_data":T[1])}
			}

		
		if(hasFeature("vep_utr") && !isBlank(params.genomes[genomeId],"ucsc_name") && ( params.genomes[genomeId].ucsc_name.equals("hg19") || params.genomes[genomeId].ucsc_name.equals("hg38"))) {
			veputr_ch = UTR_ANNOTATOR_DOWNLOAD_01([:], genomeId)
        		version_ch = version_ch.mix(veputr_ch.version)

			rows = rows.combine(veputr_ch.output).
					map{T->T[0].plus("vep_utr":T[1])}

			}


		if(hasFeature("phyloP") &&  !isBlank(params.genomes[genomeId],"phyloP_bigwig_url")) {
			phyloP_ch =  WGET_PHYLOP([:], [url:params.genomes[genomeId].phyloP_bigwig_url, output:"phyloP.bw"])

			rows = rows.combine(phyloP_ch.output).
					map{T->T[0].plus("phyloP_bigwig":T[1])}
			}

		if(hasFeature("phastCons") &&  !isBlank(params.genomes[genomeId],"phastCons_bigwig_url")) {
			phastCons_ch =  WGET_PHASTCONS([:], [url:params.genomes[genomeId].phastCons_bigwig_url, output:"phastCons.bw"])

			rows = rows.combine(phastCons_ch.output).
					map{T->T[0].plus("phastCons_bigwig":T[1])}
			}


		annot_vcf_ch = APPLY_ANNOTATION(
			[:],
			genomeId,
			rows
			)


		rows = annot_vcf_ch.output.map{T->T[0].plus("annot_vcf":T[1])}
		version_ch = version_ch.mix(annot_vcf_ch.version)

		version_ch = MERGE_VERSION("VCF annotation", version_ch.collect())
	emit:
		version= version_ch
		output = rows
	}


process complementCaptures {
memory "3g"
input:
	val(meta)
output:
	path("norm.captures.tsv"),emit:output
script:
if(isBlank(meta,captures))
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
afterScript "rm -f org.tsv genes.tsv snpeff.errors"
memory "10g"
input:
        val(meta)
        val(genomeId)
output:
	path("snpEff.config"),emit:output
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
tag "${row.vcf} ${row.interval?:row.bed?:""}"
afterScript "rm -rf TMP"
memory '5 g'
maxForks 30

input:
	val(meta)
	val(genomeId)
	val(row)
output:
	tuple val(row),path("contig.bcf"),emit:output
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


	#
	# CUSTOM VCFFILTERJDK expression
	#
	if ${!params.annotations.vcffilterjdk_args1.trim().isEmpty()} ; then

		java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar  \${JVARKIT_DIST}/jvarkit.jar vcffilterjdk \
			${params.annotations.vcffilterjdk_args1} \
			TMP/jeter1.vcf > TMP/jeter2.vcf
		
		mv -v TMP/jeter2.vcf TMP/jeter1.vcf

	fi

	# annotate variant ID
	if ${hasFeature("dbsnp_id") && genome.containsKey("dbsnp")} ; then

		# file must be indexed
		bcftools view -O b -o TMP/jeter1.bcf TMP/jeter1.vcf
		bcftools index TMP/jeter1.bcf
		
		bcftools annotate --annotations "${genome.dbsnp}"  --regions-file "TMP/jeter.123.bed"  -c ID -O v -o TMP/jeter1.vcf TMP/jeter1.bcf
		rm TMP/jeter1.bcf TMP/jeter1.bcf.csi
	fi


	#
	# MICRO ORF
	#
	if ${hasFeature("gff3_uorf") && !isBlank(genome,"gtf")} then

		java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar  \${JVARKIT_DIST}/jvarkit.jar vcfscanupstreamorf \
			--canonical --exclude-cds --gtf "${genome.gtf}" -R "${reference}" --strong TMP/jeter1.vcf > TMP/jeter2.vcf
		
		mv -v TMP/jeter2.vcf TMP/jeter1.vcf
	
	fi


	# bcftools CSQ
        if ${hasFeature("bcftools_csq") && !isBlank(genome,"gff3")}  ; then

  		bcftools csq -O v --force --local-csq --ncsq 10000 --fasta-ref "${reference}" --gff-annot "${genome.gff3}" TMP/jeter1.vcf > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi

	# annotations with bcftools annotate
	if ${row.annotate && !row.annotate.isEmpty()} ; then
		
		bcftools view -O u -o TMP/jeter1.bcf TMP/jeter1.vcf

		${row.annotate.collect{T->"####\n# ${T.name}\n####\nbcftools index -f TMP/jeter1.bcf && bcftools annotate --annotations '${T.tabix}' --columns '${T.columns}' ${T.tabix.name.endsWith(".vcf.gz") || T.tabix.name.endsWith(".bcf")?"":"--header-lines '${T.header}'"} ${T.merge_logic.isEmpty()?"":"--merge-logic '${T.merge_logic}'"}  ${T.minoverlap.equals(".")?"":"--min-overlap '${T.minoverlap}'"} -O b -o TMP/jeter2.bcf TMP/jeter1.bcf\nmv -v TMP/jeter2.bcf TMP/jeter1.bcf"}.join("\n")}

		bcftools view -O v -o TMP/jeter1.vcf TMP/jeter1.bcf
		rm -f TMP/jeter1.bcf TMP/jeter1.bcf.csi
	fi


	# polyx
	if ${hasFeature("polyx") && (params.annotations.polyx as int)>0} ; then
		java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfpolyx \
			-n '${params.annotations.polyx}' \
			--reference '${reference}' \
			--tag POLYX TMP/jeter1.vcf > TMP/jeter2.vcf

		mv -v TMP/jeter2.vcf TMP/jeter1.vcf
	fi
	

	# allelic ratio
	if ${hasFeature("allelic_ratio") && (params.annotations.allelic_ratio as double as double)>=0} ; then
		java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar  \${JVARKIT_DIST}/jvarkit.jar vcffilterjdk \
			--nocode \
			${isSoftFilter("allelic_ratio")?"--filter HET_BAD_AD_RATIO ":""} \
			-e 'return variant.getGenotypes().stream().filter(G->G.isHet() && G.hasAD() && G.getAD().length==2).allMatch(G->{int array[]=G.getAD();double r= array[1]/(double)(array[0]+array[1]);final double x=${params.annotations.allelic_ratio as double}; return (r>=x && r<=(1.0 - x)) ;});' TMP/jeter1.vcf > TMP/jeter2.vcf
		
		mv -v TMP/jeter2.vcf TMP/jeter1.vcf
	fi


	# LOW DP
	if ${hasFeature("DP") && (params.annotations.low_DP as int)>=0} ; then
		java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar   \${JVARKIT_DIST}/jvarkit.jar vcffilterjdk \
			--nocode \
			${isSoftFilter("DP")?"--filter LOW_DEPTH_params.annotations.low_DP ":""} \
			-e 'return variant.getGenotypes().stream().filter(G->G.hasDP() && G.getAlleles().stream().anyMatch(A->!(A.isNoCall() || A.isReference())) ).mapToInt(G->G.getDP()).average().orElse(1.0+${params.annotations.low_DP}) > ${params.annotations.low_DP} ;' TMP/jeter1.vcf > TMP/jeter2.vcf
	
		mv -v TMP/jeter2.vcf TMP/jeter1.vcf
	fi

	# HIGH DP
	if ${hasFeature("DP") && (params.annotations.high_DP as int)>=0} ; then
		java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar   \${JVARKIT_DIST}/jvarkit.jar vcffilterjdk \
			--nocode \
			${isSoftFilter("DP")?"--filter HIGH_DEPTH_${params.annotations.high_DP} ":""} \
			-e 'return variant.getGenotypes().stream().filter(G->G.hasDP() && G.getAlleles().stream().anyMatch(A->!(A.isNoCall() || A.isReference())) ).mapToInt(G->G.getDP()).average().orElse(-1.0) < ${params.annotations.high_DP} ;' TMP/jeter1.vcf > TMP/jeter2.vcf
	
		mv -v TMP/jeter2.vcf TMP/jeter1.vcf
	fi

	#
	# GQ
	#
	if ${hasFeature("GQ") && (params.annotations.low_GQ as int ) >=0 } ; then
		java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffilterjdk \
			--nocode \
			${isSoftFilter("GQ")?"--filter LOW_GQ_${params.annotations.low_GQ} ":""} \
			-e 'return variant.getGenotypes().stream().filter(G->G.hasGQ() && G.getAlleles().stream().anyMatch(A->!(A.isNoCall() || A.isReference())) ).allMatch(G->G.getGQ()>=${params.annotations.low_GQ});' TMP/jeter1.vcf > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf

	fi


	#
	# SNPEFF
	#
	if ${hasFeature("snpeff")} ; then
		if ${hasFeature("native_snpeff") && !isBlank(genome,"snpeff_database_name") } ; then
			module load snpEff
			java -Djava.io.tmpdir=TMP -jar "\${SNPEFF_JAR}" eff -config "\${SNPEFF_CONFIG}" \
				-nodownload -noNextProt -noMotif -noInteraction -noLog -noStats -chr chr -i vcf -o vcf '${genome.snpeff_database_name}' < TMP/jeter1.vcf > TMP/jeter2.vcf
			mv -v TMP/jeter2.vcf TMP/jeter1.vcf
		elif ${row.containsKey("snpeff_data")} ; then
			module load snpEff
			java -Djava.io.tmpdir=TMP -jar "\${SNPEFF_JAR}" eff -config '${row.snpeff_data}' \
				-nodownload -noNextProt -noMotif -noInteraction -noLog -noStats -chr chr -i vcf -o vcf '${genomeId}' < TMP/jeter1.vcf > TMP/jeter2.vcf
			mv -v TMP/jeter2.vcf TMP/jeter1.vcf
		fi
	fi

	#
	# VEP
	#
        if ${hasFeature("vep") && genome.containsKey("vep_module") && !isBlank(genome,"vep_invocation") }  ; then
		module load ${genome.vep_module}

		${row.containsKey("vep_utr")?
			"set +u; export PERL5LIB=\${PERL5LIB}:${row.vep_utr.getParent()}":
			""
		}

		${!isBlank(row,"vep_utr") && !isBlank(genome,"vep_invocation") ?
			genome.vep_invocation.replace("--cache", "--plugin UTRannotator,${row.vep_utr.toRealPath()} --cache") :
			genome.vep_invocation
			} < TMP/jeter1.vcf > TMP/jeter2.vcf

		mv -v TMP/jeter2.vcf TMP/jeter1.vcf 1>&2
		module unload ${genome.vep_module}
	fi


	####
	#
	# phastCons
	#
	#
	if ${hasFeature("phastCons") && ${!isBlank(row,"phastCons_bigwig")} : then
		java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfbigwig  \
			--bigwig '${row.phastCons_bigwig?:"NO_FILE"}' \
			TMP/jeter1.vcf > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf

	fi

	####
	#
	# phyloP
	#
	#
	if ${hasFeature("phyloP") && ${!isBlank(row,"phyloP_bigwig")} : then

		java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfbigwig  \
			--bigwig '${row.phyloP_bigwig?:"NO_FILE"}' \
			TMP/jeter1.vcf > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf

	fi



	#
	# VCFFILTERSO
	#
        if ${hasFeature("vcffilterso") && !isBlank(params.annotations,"soacn")} && bcftools view --header-only TMP/jeter1.vcf | grep -m 1 -E '^##INFO=<ID=(BCSQ|CSQ|ANN),'  ; then
		java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffilterso \
			${isSoftFilter("vcffilterso")?"--filterout  BAD_SO":""} \
			--acn "${params.annotations.soacn}" TMP/jeter1.vcf > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf

	fi


	#
	# GNOMAD GENOME
	#
	if ${hasFeature("gnomad_genome") && !isBlank(genome,"gnomad_genome")} ; then

		java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfgnomad \
			--bufferSize 10000 \
			--max-af ${params.annotations.gnomad_max_AF} \
			--gnomad "${genome.gnomad_genome}" --fields "${params.annotations.gnomad_population}" TMP/jeter1.vcf > TMP/jeter2.vcf
                mv TMP/jeter2.vcf TMP/jeter1.vcf


		if  ${!isSoftFilter("gnomad")} ; then
	
			bcftools view  -e 'FILTER ~ "GNOMAD_GENOME_BAD_AF" || FILTER ~ "GNOMAD_GENOME_RF"'  TMP/jeter1.vcf > TMP/jeter2.vcf
	                mv TMP/jeter2.vcf TMP/jeter1.vcf			

		fi
	fi

	#
	# GNOMAD EXOME
	#
	if ${hasFeature("gnomad_exome") && !isBlank(genome,"gnomad_exome")} ; then

		java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfgnomad \
			--bufferSize 10000 \
			--max-af ${params.annotations.gnomad_max_AF} \
			--gnomad "${genome.gnomad_exome}" --fields "${params.annotations.gnomad_population}" TMP/jeter1.vcf > TMP/jeter2.vcf
                mv TMP/jeter2.vcf TMP/jeter1.vcf


		if  ${!isSoftFilter("gnomad")} ; then
	
			bcftools view  -e 'FILTER ~ "GNOMAD_EXOME_BAD_AF" || FILTER ~ "GNOMAD_EXOME_RF"'  TMP/jeter1.vcf > TMP/jeter2.vcf
	                mv TMP/jeter2.vcf TMP/jeter1.vcf			

		fi

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

/*


*/ 		


process TODOO {

script:
"""

	# select with jvarkit
	if [ ! -z "${isBlank(meta.select)?"":"Y"}" ] ; then
		java  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/vcffilterjdk.jar --nocode  ${meta.select} TMP/jeter1.vcf > TMP/jeter2.vcf
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

