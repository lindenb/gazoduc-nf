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
def gazoduc = gazoduc.Gazoduc.getInstance(params).putGnomad()

include {getModules;moduleLoad} from '../../modules/utils/functions.nf'
include {WGSELECT_EXCLUDE_BED_01 } from './wgselect.exclude.bed.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_CONCAT_PER_CONTIG_01} from '../bcftools/bcftools.concat.contigs.01.nf' 



workflow WGSELECT_01 {
	take:
		meta
		genomeId
		rows 
	main:
		version_ch = Channel.empty()
		
		
		exclude_ch = WGSELECT_EXCLUDE_BED_01(genomeId)
		version_ch = version_ch.mix(exclude_ch.version)


		annotate_ch = ANNOTATE(meta, genomeId, rows.combine(exclude_ch.bed).map{T->T[0].plus(blacklisted:T[1])})

		version_ch = version_ch.mix(annotate_ch.version.first())
		
		cat_files_ch = COLLECT_TO_FILE_01([:], annotate_ch.bed_vcf.map{T->T[1]}.collect())

		concat_ch = BCFTOOLS_CONCAT_PER_CONTIG_01([suffix:".bcf"],  cat_files_ch.output)
		version_ch = version_ch.mix(concat_ch.version)



		digest_ch = DIGEST_VARIANT_LIST(annotate_ch.variants_list.collect())
		version_ch = version_ch.mix(digest_ch.first())


		version_ch = MERGE_VERSION("wgselect", version_ch.collect())
	emit:
		version = version_ch /** version */
		variants_list = digest_ch.output /** file containing count of variants at each step */
		contig_vcfs = concat_ch.vcfs /** path to all vcf concatenated per contigs */
		bed = concat_ch.bed /** path to all VCFs as bed */
		vcfs = cat_files_ch.output /** path to all chunks of vcf */
	}



process ANNOTATE {
tag "${file(row.vcf).name} ${row.bed?:row.interval}"
cache "lenient"
cpus {task.attempt}
errorStrategy 'retry'
maxRetries 5
memory "5g"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(genomeId)
	val(row)
output:
	tuple path("contig.bed"),path("contig.bcf"), emit: bed_vcf
	path("variant_list.txt.gz"), emit:variants_list
	path("version.xml"),emit:version
script:
	if(!row.containsKey("pedigree")) throw new IllegalArgumentException("row.pedigree is missing");
	if(!row.containsKey("vcf")) throw new IllegalArgumentException("row.vcf is missing");
	if(!(row.containsKey("bed") || row.containsKey("interval"))) throw new IllegalArgumentException("row.bed/interval is missing");
	if(!row.containsKey("hard_filters")) throw new IllegalArgumentException("row.hard_filters is missing");
	//
	def vcf = row.vcf
	def jvarkitped = row.pedigree
	def genome = params.genomes[genomeId]
	def reference = genome.fasta	
	def blacklisted = row.blacklisted
	def gnomadgenome = genome.gnomad_genome
	def gnomadgenomefilterexpr = genome.ucsc_name.equals("hg19")? params.wgselect.gnomadgenomefilterexpr_hg19 : ( genome.ucsc_name.equals("hg38")? params.wgselect.gnomadgenomefilterexpr_hg38 : "")
	def mapability= genome.mapability_bigwig
	def snpeffDb = genome.snpeff_database_name
	def vep_module = genome.vep_module?:""
	def vep_invocation = genome.vep_invocation?:""
	def max_alleles_count = (params.wgselect.max_alleles_count as int)
	def polyx = (params.wgselect.polyx as int)
	def gnomadPop = (params.wgselect.gnomadPop)
	def gnomadAF = (params.wgselect.gnomadAF as double)
	def soacn = params.wgselect.soacn
	def exclude_soacn = params.wgselect.exclude_soacn?:""
	def inverse_so = (params.wgselect.inverse_so as boolean)
	def f_missing = (params.wgselect.f_missing as double)
	def maxmaf = (params.wgselect.maxmaf as double)
	def fisherh = (params.wgselect.fisherh as double)
	def minDP= (params.wgselect.minDP as int)
	def maxDP= (params.wgselect.maxDP as int)
	def lowGQ =( params.wgselect.lowGQ as int)
	def hwe = (params.wgselect.hwe as double)
	def minGQsingleton = (params.wgselect.minGQsingleton as int)
	def minRatioSingleton  = (params.wgselect.minRatioSingleton as double)
	def annot_method = (params.wgselect.annot_method)
	def cadd_phred = (params.wgselect.cadd_phred as double)
	def cadd_tabix = genome.cadd_tabix
	def apply_hard_filters_arguments = row.hard_filters

"""
hostname 1>&2
${moduleLoad("jvarkit bcftools bedtools")}
set -x
mkdir -p TMP
touch TMP/variant_list.txt
which java 1>&2
which javac 1>&2
javac -version 1>&2
echo "\${JAVA_HOME}"


	if ${row.containsKey("bed")} ; then
		ln -s '${row.bed}' TMP/contig.bed
	elif ${row.containsKey("interval")} ; then
		echo '${row.interval}' | awk -F '[:-]' '{printf("%s\t%s\t%s\\n",\$1,int(\$2)-1,\$3);}' > TMP/contig.bed
	else
		echo "No interval or bed defined" 1>&2
		exit -1
	fi


	function countIt {
		if ${params.wgselect.with_count as boolean} ; then
			echo -n "\$1\t" >> TMP/variant_list.txt
			bcftools query -f '%CHROM:%POS:%REF:%ALT\\n' "\$2" | sed 's/^chr//' | LC_ALL=C sort -T . | uniq  > TMP/tmp.A.txt
			bcftools query -f '%CHROM:%POS:%REF:%ALT\\n' "\$3" | sed 's/^chr//' | LC_ALL=C sort -T . | uniq  > TMP/tmp.B.txt
			LC_ALL=C comm -23 TMP/tmp.A.txt TMP/tmp.B.txt > TMP/tmp.C.txt
			cat TMP/tmp.A.txt | wc -l | tr "\\n" "\t" >> TMP/variant_list.txt
			cat TMP/tmp.B.txt | wc -l | tr "\\n" "\t" >> TMP/variant_list.txt
			cat TMP/tmp.C.txt | wc -l | tr "\\n" "\t" >> TMP/variant_list.txt
			cat TMP/tmp.C.txt | paste -d ';' -s >> TMP/variant_list.txt
			rm TMP/tmp.A.txt TMP/tmp.B.txt TMP/tmp.C.txt
		fi
	}

	cat <<-EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">wgselect for one bed</entry>
		<entry key="vcf">${vcf}</entry>
		<entry key="ped">${jvarkitped}</entry>
		<entry key="steps">
	EOF

	## Extract case/controls from pedigree
	if [ "${jvarkitped.name}" == "NO_FILE" ] ; then
		touch TMP/cases.txt TMP/controls.txt

		cat <<-EOF >> version.xml
		<properties>
			<entry key="description">extract case/controls from pedigree. No pedigree was provided, so there is none.</entry>
		</properties>
		EOF
		
	else
		# extract list of samples cases and controls
		awk -F '\t' '(\$6=="case" ||  \$6=="affected") {print \$2;}' "${jvarkitped}" | sort | uniq > TMP/cases.txt
		awk -F '\t' '(\$6=="control" ||  \$6=="unaffected") {print \$2;}' "${jvarkitped}" | sort | uniq > TMP/controls.txt
		cat <<-EOF >> version.xml
		<properties>
			<entry key="description">extract case/controls from pedigree</entry>
			<entry key="cases.count">\$(wc -l < TMP/cases.txt)</entry>
			<entry key="controls.count">\$(wc -l < TMP/controls.txt)</entry>
		</properties>
		EOF
	fi

	# blacklisted region overlapping #####################################################################"
	bedtools intersect -a "TMP/contig.bed" -b "${blacklisted}" > TMP/jeter.blacklisted.bed
	# prevent empty file
	if [ ! -s TMP/jeter.blacklisted.bed ] ; then
		tail -1 "${reference}.fai" | awk '{printf("%s\t0\t1\\n",\$1);}' > TMP/jeter.blacklisted.bed
	fi
	cat <<-EOF >> version.xml
	<properties>
		<entry key="description">excluded regions are reduced by overlap over bed file using bedtools</entry>
		<entry key="vcf">${vcf}</entry>
	</properties>
	EOF
	
	

	# extract variants ######################################################################################
	if ${vcf.endsWith(".list")} ; then
		bcftools concat --file-list "${vcf}" --regions-file "TMP/contig.bed" -O u  --allow-overlaps --remove-duplicates -o TMP/jeter1.bcf
	else
		bcftools view  --regions-file "TMP/contig.bed" -O u -o TMP/jeter1.bcf "${vcf}"
	fi	

	cat <<- EOF >> version.xml
	<properties>
		<entry key="description">extract variants from main vcf in bed region</entry>
		<entry key="vcf">${vcf}</entry>
	</properties>
	EOF


	if ${params.wgselect.with_count as boolean} ; then
		bcftools query -f '.\\n' TMP/jeter1.bcf | awk 'END{printf("initial\t%s\t0\t0\t\\n",NR);}' >> TMP/variant_list.txt
	fi



	# gatk hard filtering #############################################################
	if ${!file(apply_hard_filters_arguments).name.equals("NO_FILE") } ; then
		## conflic betwwen java for jvarkit and gatk "Duplicate cpuset controllers detected" on stdout
		module unload jvarkit
		module load gatk/0.0.0
		bcftools view -O z -o TMP/jeter1.vcf.gz TMP/jeter1.bcf
		bcftools index -f -t TMP/jeter1.vcf.gz

		gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" VariantFiltration \
	        	-L "TMP/contig.bed" \
		        -V 'TMP/jeter1.vcf.gz' \
	        	-R '${reference}' \
		        -O TMP/jeter2.vcf.gz \
        		--arguments_file ${apply_hard_filters_arguments}

		bcftools view --apply-filters 'PASS,.'  -O b -o TMP/jeter1.bcf TMP/jeter2.vcf.gz
		countIt "gatkHardFilters" TMP/jeter1.vcf.gz TMP/jeter2.vcf.gz

		rm TMP/jeter2.vcf.gz TMP/jeter2.vcf.gz.tbi TMP/jeter1.vcf.gz TMP/jeter1.vcf.gz.tbi

		module unload gatk/0.0.0
		module load jvarkit


		cat <<- EOF >> version.xml
		<properties>
			<entry key="description">remove variant failing gatk hard filters</entry>
		</properties>
		EOF
		
	fi


	# min max alleles  ############################################################################
	bcftools view  -m2 -M ${max_alleles_count} -O u -o TMP/jeter2.bcf  TMP/jeter1.bcf
	countIt "too_many_alts" TMP/jeter1.bcf TMP/jeter2.bcf
	mv TMP/jeter2.bcf TMP/jeter1.bcf

	cat <<- EOF >> version.xml
	<properties>
		<entry key="description">select variants on min/max number of alleles (diallelic is 2)</entry>
		<entry key="min-alleles">2 (including REF)</entry>
		<entry key="max-alleles">${max_alleles_count} (use <code>--max_alleles_count 'x'</code> to change this)</entry>
	</properties>
	EOF


	# remove in blaclisted regions ############################################################################
	bcftools view  --targets-overlap 2 --targets-file ^TMP/jeter.blacklisted.bed -O u -o TMP/jeter2.bcf TMP/jeter1.bcf
	countIt "blackListedRegions" TMP/jeter1.bcf TMP/jeter2.bcf
	mv TMP/jeter2.bcf TMP/jeter1.bcf
	
	cat <<- EOF >> version.xml
	<properties>
		<entry key="description">Ignore variants overlapping blacklisted region</entry>
		<entry key="exclude bed">${blacklisted}</entry>
	</properties>
	EOF


	# remove all annotations ################################################################################
	bcftools annotate --force -x '^INFO/AC,INFO/AN,INFO/ReadPosRankSum,INFO/MQRankSum,INFO/MQ,INFO/FS,INFO/QD,INFO/SOR,FILTER' -O u -o TMP/jeter2.bcf TMP/jeter1.bcf
	countIt "annotateX" TMP/jeter1.bcf TMP/jeter2.bcf
	mv TMP/jeter2.bcf TMP/jeter1.bcf


	cat <<- EOF >> version.xml
	<properties>
		<entry key="description">remove annotations from original VCF</entry>
		<entry key="bcftools.annotate.x">^INFO/AC,INFO/AN,INFO/ReadPosRankSum,INFO/MQRankSum,INFO/MQ,INFO/FS,INFO/QD,INFO/SOR,FILTER</entry>
	</properties>
	EOF


	# could happen that some variant are discarded here: saw some gatk4 variants where *NO* genotype was called. #############################
	if test -s TMP/cases.txt || test -s TMP/controls.txt ; then
		cat TMP/cases.txt TMP/controls.txt | sort -T TMP | uniq > TMP/jeter.samples
		bcftools view --trim-alt-alleles --samples-file 'TMP/jeter.samples' -O u -o TMP/jeter2.bcf TMP/jeter1.bcf
		rm TMP/jeter.samples
		countIt "samples" TMP/jeter1.bcf TMP/jeter2.bcf
		mv TMP/jeter2.bcf TMP/jeter1.bcf

		cat <<- EOF >> version.xml
		<properties>
			<entry key="description">remove unused samples. Remove unused alleles. Saw some gatk4 variants where *NO* genotype was called</entry>
		</properties>
		EOF
	fi


	# force recalculation of AF/AC/AN #########################################################
	bcftools  +fill-tags -O u -o TMP/jeter2.bcf TMP/jeter1.bcf -- -t AN,AC,AF
	countIt "filltags" TMP/jeter1.bcf TMP/jeter2.bcf
	mv TMP/jeter2.bcf TMP/jeter1.bcf

	cat <<- EOF >> version.xml
	<properties>
		<entry key="description">fill tags in INFO using <code>bcftools  +fill-tags</code></entry>
		<entry key="tags">AC,AN,AF</entry>
	</properties>
	EOF



	## too many no-call ("genotyping rate" >= 95%) #################################################
	bcftools view -i 'CHROM=="Y" || CHROM=="chrY" ||  F_MISSING < ${f_missing}' -O b -o TMP/jeter2.bcf TMP/jeter1.bcf
	countIt "F_MISSING_${f_missing}" TMP/jeter1.bcf TMP/jeter2.bcf

	cat <<- EOF >> version.xml
	<properties>
		<entry key="description">remove to many NO_CALL <code>./.</code></entry>
		<entry key="when">CHROM!=Y or F_MISSING &lt;= ${f_missing})</entry>
	</properties>
	EOF

	mv TMP/jeter2.bcf TMP/jeter1.bcf
	
	## sex et homvar (1 homvar and 0 het)
	cat <<- EOF >> version.xml
	<properties>
                <entry key="description">remove variant on autosome if no HET and found at least one HOM_VAR.</entry>
	EOF

	if ${params.wgselect.with_homvar as boolean} ; then
		bcftools view -O v -o TMP/jeter2.vcf TMP/jeter1.bcf
		java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfpar TMP/jeter2.vcf > TMP/jeter1.vcf
		bcftools view -e 'INFO/SEX=0 &&  COUNT(GT="RA")==0 && COUNT(GT="AA")>0' -O u -o TMP/jeter2.bcf TMP/jeter1.vcf
		countIt "homvar" TMP/jeter1.bcf TMP/jeter2.bcf
		mv TMP/jeter2.bcf TMP/jeter1.bcf
		rm TMP/jeter1.vcf TMP/jeter2.vcf

	cat <<- EOF >> version.xml
		<entry key="enabled">true</entry>
	</properties>
	EOF

	else

	cat <<- EOF >> version.xml
		<entry key="enabled">false</entry>
	</properties>
	EOF

	fi

	# split multiallelic #################################################################################""
	cat <<- EOF >> version.xml
	<properties>
                <entry key="description">Split multiallelic with <code>bcftools norm --multiallelics --both</code></entry>
	EOF


	if [ "${max_alleles_count}" != "2" ] ; then

		bcftools norm -f "${reference}" --multiallelics -both  -O u -o TMP/jeter2.bcf TMP/jeter1.bcf
		countIt "norm" TMP/jeter1.bcf TMP/jeter2.bcf
		mv TMP/jeter2.bcf TMP/jeter1.bcf

		echo '<entry key="enabled">true</entry></properties>' >> version.xml

	else

		echo '<entry key="enabled">false</entry></properties>' >> version.xml
	fi

	# not in pedigree ###########################################################################
	bcftools view -i 'AC[*]>0' -O v -o TMP/jeter2.vcf TMP/jeter1.bcf 

	countIt "AC_GT_0" TMP/jeter1.bcf TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf
	rm TMP/jeter1.bcf

	cat <<- EOF >> version.xml
	<properties>
		<entry key="name">remove-samples-not-in-pedigree</entry>
		<entry key="description">after some samples have been removed (e.g: there in the vcf but not in cases/ctrls). There can be some variants without any ALT allele because that will be removed here.</entry>
	</properties>
	EOF


	# ignore spanning deletions #################################################################
	awk -F '\t' '(\$0 ~ /^#/ || \$5!="*")'  TMP/jeter1.vcf > TMP/jeter2.vcf


	cat <<- EOF >> version.xml
	<properties>
		<entry key="name">ignore spanning deletions.</entry>
		<entry key="description">after <code>bcftools norm</code>, we remove variants where the only allele is <code>&lt;*&gt;</code>. See <a>https://gatk.broadinstitute.org/hc/en-us/articles/360035531912-Spanning-or-overlapping-deletions-allele-</a></entry>
	</properties>
	EOF



	countIt "spandel" TMP/jeter1.vcf TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf


	## polyx ###################################################################################

	cat <<- EOF >> version.xml
	<properties>
		<entry key="description">remove variant near a poly-x</entry>
		<entry key="poly-x">${polyx}</entry>
	EOF

	if ${ (polyx as int) > 1 } ; then
		echo "##show me vcf" 1>&2
		cat TMP/jeter1.vcf | head 1>&2

		java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar vcfpolyx -R "${reference}" --tag POLYX -n "${polyx}" TMP/jeter1.vcf |\
		bcftools view -e 'FILTER~"POLYX_ge_${polyx}"' > TMP/jeter2.vcf
		countIt "polyx${polyx}" TMP/jeter1.vcf TMP/jeter2.vcf

		echo '<entry key="enabled">true</entry></properties>' >> version.xml

		mv TMP/jeter2.vcf TMP/jeter1.vcf
	else

		echo '<entry key="enabled">false</entry></properties>' >> version.xml
	fi
	
	## CADD ######################################################################################

	cat <<- EOF >> version.xml
	<properties>
		<entry key="description">remove variant with low CADD phred score</entry>
	EOF

	if ${!cadd_tabix.isEmpty() && (cadd_phred as double) > 0}  ; then
        	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar vcfcadd \
			--tabix "${cadd_tabix}" TMP/jeter1.vcf > TMP/jeter2.vcf
	      	mv TMP/jeter2.vcf TMP/jeter1.vcf
		
		bcftools view -e 'INFO/CADD_PHRED < ${cadd_phred as double}' -O v -o TMP/jeter2.vcf TMP/jeter1.vcf

		countIt "CADD" TMP/jeter1.vcf TMP/jeter2.vcf
	      	mv TMP/jeter2.vcf TMP/jeter1.vcf

		cat <<- EOF >> version.xml
			<entry key="cadd.phred">${cadd_phred}</entry>
			<entry key="cadd.file">${cadd_tabix}</entry>
			<entry key="enabled">true</entry>
		</properties>
		EOF


	else

		echo '<entry key="enabled">false</entry></properties>' >> version.xml

	fi





	# test MAF

	cat <<- EOF >> version.xml 
	<properties>
		<entry key="name">internal maf</entry>
		<entry key="description">remove variant if internal MAF is too high. Nextflow parameter is <code>--maxmaf (value [0.0 -  1.0] )</code> </entry>
		<entry key="max MAF">${maxmaf}</entry>
	EOF
	
	if ${maxmaf>=0} && test -s TMP/cases.txt && test -s TMP/controls.txt ; then
		java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfburdenmaf \
			--pedigree "${jvarkitped}" --prefix "" --min-maf 0  --max-maf "${maxmaf}"  TMP/jeter1.vcf   > TMP/jeter2.vcf
		countIt "MAF" TMP/jeter1.vcf TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf

		echo '<entry key="enabled">true</entry></properties>' >> version.xml
	else
	
		echo '<entry key="enabled">false (or no case/controls)</entry></properties>' >> version.xml

	fi

	# fisher per variant
	cat <<- EOF >> version.xml
        <properties>
                <entry key="name">fisher H</entry>
		<entry key="description">remove variant if fisher test per variant is tool low</entry>
		<entry key="min-fisher">${fisherh}</entry>
		<entry key="max-fisher">1.0</entry>
	EOF
	
	if ${fisherh >= 0.0} && test -s TMP/cases.txt && test -s TMP/controls.txt ; then
		java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfburdenfisherh --filter '' --pedigree "${jvarkitped}" --min-fisher "${fisherh}"  TMP/jeter1.vcf   > TMP/jeter2.vcf
		countIt "fisherH" TMP/jeter1.vcf TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf
		bcftools annotate -x 'FILTER/CTRL_CASE_RATIO' TMP/jeter1.vcf > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf

		echo '<entry key="enabled">true</entry></properties>' >> version.xml

	else

		echo '<entry key="enabled">false (or not case/control)</entry></properties>' >> version.xml

	fi

	# low or high DP

	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffilterjdk --nocode  -e 'final double dp= variant.getGenotypes().stream().filter(G->G.isCalled() && G.hasDP()).mapToInt(G->G.getDP()).average().orElse(${minDP}); return dp>=${minDP} && dp<=${maxDP};'  TMP/jeter1.vcf > TMP/jeter2.vcf
	countIt "lowDP" TMP/jeter1.vcf TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf


	cat <<- EOF >> version.xml
        <properties>
                <entry key="name">depth</entry>
		<entry key="description">remove variant if INFO/DP per variant is tool low / too high</entry>
		<entry key="min DP">${minDP}</entry>
		<entry key="max DP">${maxDP}</entry>
	</properties>
	EOF


	# all genotypes with ALT must have GQ >= 'x'
	cat <<- EOF >> version.xml
        <properties>
                <entry key="name">LOW GQ</entry>
		<entry key="description">ALL genotypes carrying a ALT must be GQ >= ${lowGQ}</entry>
		<entry key="GQ">${lowGQ}</entry>
	EOF


	if [ "${lowGQ}" -gt 0 ] ; then
		## low GQ all genotypes carrying a ALT must be GQ > 'x'
		java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffilterjdk --nocode  -e 'return variant.getGenotypes().stream().filter(g->g.isCalled() && !g.isHomRef() && g.hasGQ()).allMatch(g->g.getGQ()>=${lowGQ});' TMP/jeter1.vcf > TMP/jeter2.vcf
		countIt "lowGQ${lowGQ}" TMP/jeter1.vcf TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf

		echo '<entry key="enabled">true</entry></properties>' >> version.xml

	else
		
		echo '<entry key="enabled">false</entry></properties>' >> version.xml

	fi

	## singleton
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffilterjdk --nocode  -e 'Genotype singleton=null; for(final Genotype g: variant.getGenotypes()) {if(g.isCalled() && !g.isHomRef()) { if(singleton!=null) return true;singleton=g;}} if(singleton!=null && singleton.isFiltered()) return false; if(singleton!=null && singleton.isHet() && singleton.hasGQ() && singleton.getGQ()<${minGQsingleton}) return false; if(singleton !=null && singleton.hasAD() && singleton.isHet() && singleton.getAD().length==2) {int array[]=singleton.getAD();double r= array[1]/(double)(array[0]+array[1]);if(r< ${minRatioSingleton} || r>(1.0 - ${minRatioSingleton})) return false;} return true; ' TMP/jeter1.vcf > TMP/jeter2.vcf
	countIt "singleton" TMP/jeter1.vcf TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf


	cat <<- EOF >> version.xml
        <properties>
                <entry key="name">het singleton</entry>
		<entry key="description">remove variant if singleton has bad quality</entry>
		<entry key="singleton-GQ">${minGQsingleton}</entry>
		<entry key="AD-ratio-low">${minRatioSingleton}</entry>
		<entry key="AD-ratio-high">${1.0 - (minRatioSingleton as Double)}</entry>
	</properties>
	EOF



	if ${!mapability.isEmpty()} ; then

	cat <<- EOF >> version.xml
        <properties>
                <entry key="name">mapability</entry>
		<entry key="description">remove variant if singleton overlaps region of low mapability</entry>
		<entry key="file">${mapability}</entry>
		<entry key="treshold">1.0</entry>
	EOF

		# DukeMapability
		java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar vcfbigwig --aggregate avg -tag mapability \
			-B "${mapability}" TMP/jeter1.vcf > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf

		java  -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffilterjdk --nocode   \
				-e 'final String tag= "mapability"; return !variant.hasAttribute(tag) ||  (variant.getAttributeAsDouble(tag,0.0)==1.0);' TMP/jeter1.vcf > TMP/jeter2.vcf
		countIt "mapability" TMP/jeter1.vcf TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf

		echo '<entry key="enabled">true</entry></properties>' >> version.xml
	else

		echo '<entry key="enabled">false</entry></properties>' >> version.xml

	fi

	### GNOMAD GENOME #####################################################################################

	cat <<- EOF >> version.xml
        <properties>
                <entry key="name">gnomad genome variant</entry>
		<entry key="description">remove variant if frequent in gnomad genome</entry>
		<entry key="population">${gnomadPop}</entry>
		<entry key="AF">${gnomadAF}</entry>
	EOF

	if [ ! -z "${genome.gnomad_genome}"  ] ; then

		test -f "${genome.gnomad_genome}"

        	# gnomad genome
        	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfgnomad --bufferSize 1000 -F '${gnomadPop}' \
                	-g "${genome.gnomad_genome}" --max-af "${gnomadAF}" TMP/jeter1.vcf   > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf
		bcftools view -e '${gnomadgenomefilterexpr}' -O v -o TMP/jeter2.vcf TMP/jeter1.vcf
		countIt "gnomadgenome.${gnomadPop}" TMP/jeter1.vcf TMP/jeter2.vcf
	
		mv TMP/jeter2.vcf TMP/jeter1.vcf

		echo '<entry key="enabled">true</entry></properties>' >> version.xml

	else

		echo '<entry key="enabled">false</entry></properties>' >> version.xml

	fi

	### GNOMAD EXOME #####################################################################################

	cat <<- EOF >> version.xml
        <properties>
                <entry key="name">gnomad exome variant</entry>
		<entry key="description">remove variant if frequent in gnomad exome</entry>
		<entry key="population">${gnomadPop}</entry>
		<entry key="AF">${gnomadAF}</entry>
	EOF

	if [ ! -z "${genome.gnomad_exome}"  ] ; then

		test -f "${genome.gnomad_exome}"

        	# gnomad exome
        	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfgnomad   --bufferSize 1000 -F '${gnomadPop}' \
                	-g "${genome.gnomad_exome}" --max-af "${gnomadAF}" TMP/jeter1.vcf   > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf
		bcftools view -e 'FILTER~"GNOMAD_EXOME_AC0"|| FILTER~"GNOMAD_EXOME_BAD_AF"|| FILTER~"GNOMAD_EXOME_InbreedingCoeff"|| FILTER~"GNOMAD_EXOME_RF"' -O v -o TMP/jeter2.vcf TMP/jeter1.vcf
		countIt "gnomadexome.${gnomadPop}" TMP/jeter1.vcf TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf

		echo '<entry key="enabled">true</entry></properties>' >> version.xml
	else
		
		echo '<entry key="enabled">false</entry></properties>' >> version.xml

	fi

 	## FUNCTIONNAL ANNOTATION ##############################################################################
	cat <<- EOF >> version.xml
        <properties>
                <entry key="name">Functional annotations</entry>
		<entry key="method">${annot_method}</entry>
		<entry key="Sequence ontology">${soacn}</entry>
		<entry key="Exclude Sequence ontology">${exclude_soacn}</entry>
	EOF


    if [ ! -z "${annot_method}" ] ; then

	    if [  "${annot_method.toLowerCase()}" == "vep"  ]  && [ ! -z "${vep_module}" ] ; then

		module load ${vep_module}
		${vep_invocation} --output_file TMP/jeter2.vcf < TMP/jeter1.vcf
 		mv TMP/jeter2.vcf TMP/jeter1.vcf
		module unload ${vep_module}
	
		echo '<entry key="vep version">${vep_module}</entry>' >> version.xml


	    elif  [  "${annot_method.toLowerCase()}" == "snpeff" ] && [ ! -z "${snpeffDb}" ] ; then 

	   	 # snpeff
		 ${moduleLoad("snpEff/0.0.0")}
	         java  -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar "\${SNPEFF_JAR}" eff -config "\${SNPEFF_CONFIG}" -interval "TMP/contig.bed" -nodownload -noNextProt -noMotif -noInteraction -noLog -noStats -chr chr -i vcf -o vcf "${snpeffDb}" TMP/jeter1.vcf > TMP/jeter2.vcf
	         module unload ${getModules("snpEff/0.0.0")}
	         mv TMP/jeter2.vcf TMP/jeter1.vcf

		echo '<entry key="snpeff">${snpeffDb}</entry>' >> version.xml


	    else
		echo "undefined annotation method ${annot_method}:${snpeffDb}" 1>&2
		exit -1
	    fi
	


	    if ${!soacn.isEmpty()} ; then
	    
	    	# filter prediction
		java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffilterso \
			${inverse_so?"--invert":""} \
			--remove-attribute  --rmnoatt \
			--acn "${soacn}" \
		   	TMP/jeter1.vcf   2> /dev/null > TMP/jeter2.vcf

		cat <<- EOF
		<entry key="sequence.ontology">
			<properties>
				<entry key="terms">${soacn}</entry>
			</properties>
		</entry>
		EOF
		countIt "prediction" TMP/jeter1.vcf TMP/jeter2.vcf
	        mv TMP/jeter2.vcf TMP/jeter1.vcf

	    fi


	    if ${!exclude_soacn.isEmpty()} ; then
	    
	    	# filter prediction
		java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffilterso \
			--remove-attribute  --rmnoatt \
			--invert \
			--acn "${exclude_soacn}" \
		   	TMP/jeter1.vcf   2> /dev/null > TMP/jeter2.vcf

		cat <<- EOF
		<entry key="exclude_sequence.ontology">
			<properties>
				<entry key="terms">${exclude_soacn}</entry>
			</properties>
		</entry>
		EOF
		countIt "exclude.prediction" TMP/jeter1.vcf TMP/jeter2.vcf
	        mv TMP/jeter2.vcf TMP/jeter1.vcf

	    fi


		echo "</properties>" >> version.xml

    fi


    ## Hardy-Weinberg exact test

	cat <<- EOF >> version.xml
        <properties>
                <entry key="name">Hardy-Weinberg</entry>
		<entry key="description">hardy-weinberg use variants with HW  on autosomal chromosomes</entry>
	EOF


    if ${hwe >=0} ; then
	module load ${getModules("vcftools")}
	
	#set -x 
	awk '!(\$1 ~ /^(chr)?[XY]\$/ )' "TMP/contig.bed" > TMP/autosomes.bed
	awk '(\$1 ~ /^(chr)?[XY]\$/ )' "TMP/contig.bed" > TMP/sex.bed
	
	# save variants on sexual chromosomes
	if test -s TMP/sex.bed ; then
		bcftools view --targets-file TMP/sex.bed -O b -o TMP/jeter.sex.bcf TMP/jeter1.vcf
	else
		bcftools view -O b -o TMP/jeter.sex.bcf --header-only TMP/jeter1.vcf
	fi
	bcftools index TMP/jeter.sex.bcf


	# use variants with HW  on autosomal chromosomes
	if test -s TMP/autosomes.bed ; then
		bcftools view -O v -o TMP/autosomes.vcf --targets-file TMP/autosomes.bed TMP/jeter1.vcf
		vcftools --vcf TMP/autosomes.vcf --hwe "${hwe}"  --recode-INFO-all --recode --out TMP/autosome2

		bcftools view -O b -o TMP/jeter.autosomes.bcf TMP/autosome2.recode.vcf
		rm TMP/autosome2.recode.vcf TMP/autosomes.vcf
	else
		bcftools view -O b -o TMP/jeter.autosomes.bcf --header-only TMP/jeter1.vcf
	fi

	bcftools index TMP/jeter.autosomes.bcf

	bcftools concat --allow-overlaps --remove-duplicates -O v -o TMP/jeter2.vcf TMP/jeter.autosomes.bcf TMP/jeter.sex.bcf
	countIt "HW" TMP/jeter1.vcf TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf


		echo '<entry key="enabled">true</entry></properties>' >> version.xml

    else

		echo '<entry key="enabled">false</entry></properties>' >> version.xml

    fi


    # CONTRAST #############################################################################################

	cat <<- EOF >> version.xml
        <properties>
                <entry key="name">contrast</entry>
		<entry key="description">add <code>PASSOC,FASSOC,NASSOC,NOVELAL,NOVELGT</code> with <code>bcftools +contrast</code></entry>
	EOF
	


	if ${(params.wgselect.with_contrast as boolean)}  && test -s TMP/cases.txt && test -s TMP/controls.txt ; then
		bcftools +contrast \
			-0 "TMP/controls.txt" \
			-1 "TMP/cases.txt" \
			-a PASSOC,FASSOC,NASSOC,NOVELAL,NOVELGT -O v -o TMP/jeter2.vcf TMP/jeter1.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf


		echo '<entry key="enabled">true</entry></properties>' >> version.xml

	else

		echo '<entry key="enabled">false</entry></properties>' >> version.xml

	fi	


bcftools view  -O u TMP/jeter1.vcf |\
bcftools sort -T TMP --max-mem "${task.memory.giga}G" -O b -o "contig.bcf" 
bcftools index  contig.bcf


mv TMP/contig.bed ./

echo '</entry></properties>' >> version.xml

# check XML is OK
xmllint --noout version.xml

countIt "final" TMP/jeter1.vcf contig.bcf
gzip --best TMP/variant_list.txt
mv TMP/variant_list.txt.gz ./
"""
}


process mkFileList {
executor "local"
tag "N=${L.size()}"
input:
	val(meta)
	val(L)
output:
	path("${prefix}bed.vcf.list"),emit:list
script:
	prefix = meta.prefix?:""
"""
cat << EOF > "${prefix}bed.vcf.list"
${L.join("\n")}
EOF
"""
}

process DIGEST_VARIANT_LIST {
	tag "N=${L.size()}"
	
	input:
        	val(L)
        output:
                path("${params.prefix?:""}wgselect.count.tsv"),emit:output
                path("version.xml"),emit:version
        script:
        	prefix = params.prefix?:""
        """
        hostname 1>&2
        ${moduleLoad("datamash")}

        echo "#FILTER\tIN\tOUT\tDIFF" > "${prefix}filters.count.tsv"
        cat ${L.join(" ")} | gunzip -c | cut -f 1-4 | sort -T . -t '\t' -k1,1 |\
        	datamash  -g 1  sum 2 sum 3 sum 4 >> "${prefix}wgselect.count.tsv"
        	
	#######################################################################################"
        cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Merging the count of ${L.size()} file(s)</entry>
		<entry key="Output file">${prefix}wgselect.count.tsv</entry>
	</properties>
	EOF
        """
        }


