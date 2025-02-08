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
include {getModules;moduleLoad} from '../../modules/utils/functions.nf'
include {WGSELECT_EXCLUDE_BED_01 } from './wgselect.exclude.bed.01.nf'
include {BCFTOOLS_CONCAT} from '../bcftools/concat/main'
include {SNPEFF_DOWNLOAD} from '../../modules/snpeff/download'


workflow WGSELECT_01 {
	take:
		reference
		rows
		hard_filters
		pedigree
	main:
		
		
		exclude_ch = WGSELECT_EXCLUDE_BED_01(reference)

		snpeff_sb = SNPEFF_DOWNLOAD(params.snpeff_db)

		
		rows.view()

		annotate_ch = ANNOTATE(reference, rows , snpeff_sb.output, pedigree, hard_filters, exclude_ch.bed)

		

		concat_ch = BCFTOOLS_CONCAT( annotate_ch.output.map{[it[1],it[2]]}, file("NO_FILE"))



		digest_ch = DIGEST_VARIANT_LIST(annotate_ch.variants_list.collect())


	emit:
		//variants_list = digest_ch.output /** file containing count of variants at each step */
		//contig_vcfs = concat_ch.vcfs /** path to all vcf concatenated per contigs */
		//bed = concat_ch.bed /** path to all VCFs as bed */
		// vcfs = cat_files_ch.output /** path to all chunks of vcf */
		vcfs = Channel.empty()
	}



process ANNOTATE {
label "process_low"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
tag "${vcf.name} ${interval}"
afterScript "rm -rf TMP"
input:
	path(genome)
	tuple val(interval),path(vcf),path(vcfidx)
	path(snpEffDir)
	path(pedigree)
	path(apply_hard_filters_arguments)
	path(blacklisted)
output:
	tuple val(interval),path("*contig.bcf"),path("*contig.bcf.csi"), emit: output
	path("*variant_list.txt.gz"), emit:variants_list
script:
	def reference = genome.find{it.name.endsWith("a")}
	def gnomadgenome = params.gnomad
	def mapability= params.mapability_bigwig
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
	def cadd_phred = (params.wgselect.cadd_phred as double)
	def cadd_tabix = params.cadd

"""
hostname 1>&2
set -x
mkdir -p TMP
touch TMP/variant_list.txt


# jvarkit executable in conda
JD1=`which jvarkit`
echo "JD1=\${JD1}" 1>&2
# directory of jvarkit
JD2=`dirname "\${JD1}"`
# find the jar itself
JVARKIT_JAR=`find "\${JD2}/../.." -type f -name "jvarkit.jar" | head -n1`
JVARKIT_DIST=`dirname "\${JVARKIT_JAR}"`

test ! -z "\${JVARKIT_JAR}"



	echo '${interval}' | awk -F '[:-]' '{printf("%s\t%s\t%s\\n",\$1,int(\$2)-1,\$3);}' > TMP/contig.bed


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

	## Extract case/controls from pedigree
	if ${!pedigree.name.contains(".")} ; then
		touch TMP/cases.txt TMP/controls.txt TMP/pedigree.ped
	else
		# extract list of samples cases and controls
		awk -F '\t' '(\$2=="case" ||  \$2=="affected") {print \$1;}' "${pedigree}" | sort | uniq > TMP/cases.txt
		awk -F '\t' '(\$2=="control" ||  \$2=="unaffected") {print \$1;}' "${pedigree}" | sort | uniq > TMP/controls.txt
		awk -F '\t' '{printf("%s\t%s\t0\t0\t0\tcase\\n",\$1,\$1);}' TMP/cases.txt >  TMP/pedigree.ped
		awk -F '\t' '{printf("%s\t%s\t0\t0\t0\tcontrol\\n",\$1,\$1);}' TMP/controls.txt >> TMP/pedigree.ped
	fi

	# blacklisted region overlapping #####################################################################"
	bedtools intersect -a "TMP/contig.bed" -b "${blacklisted}" -nonamecheck > TMP/jeter.blacklisted.bed
	# prevent empty file
	if [ ! -s TMP/jeter.blacklisted.bed ] ; then
		tail -1 "${reference}.fai" | awk '{printf("%s\t0\t1\\n",\$1);}' > TMP/jeter.blacklisted.bed
	fi
	
	

	# extract variants  ######################################################################################
	if ${vcf.endsWith(".list")} ; then

		bcftools concat --file-list "${vcf}" --regions-file "TMP/contig.bed" -O u  --allow-overlaps --remove-duplicates -o TMP/jeter1.bcf
		bcftools view ${params.wgselect.bcftools_options} -O u -o TMP/jeter2.bcf  TMP/jeter1.bcf
		mv -v TMP/jeter2.bcf TMP/jeter1.bcf
	else

		bcftools view ${params.wgselect.bcftools_options} --regions-file "TMP/contig.bed" -O u -o TMP/jeter1.bcf "${vcf}"

	fi	

	if ${params.wgselect.with_count as boolean} ; then
		bcftools query -f '.\\n' TMP/jeter1.bcf | awk 'END{printf("initial\t%s\t0\t0\t\\n",NR);}' >> TMP/variant_list.txt
	fi



	# gatk hard filtering #############################################################
	if ${apply_hard_filters_arguments.name.contains(".")} ; then
		## conflic betwwen java for jvarkit and gatk "Duplicate cpuset controllers detected" on stdout
		bcftools view -O z -o TMP/jeter1.vcf.gz TMP/jeter1.bcf
		bcftools index -f -t TMP/jeter1.vcf.gz

		gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" VariantFiltration \\
	        	-L "TMP/contig.bed" \\
		        -V 'TMP/jeter1.vcf.gz' \\
	        	-R '${reference}' \\
		        -O TMP/jeter2.vcf.gz \\
        		--arguments_file "${apply_hard_filters_arguments}"

		bcftools view --apply-filters 'PASS,.'  -O b -o TMP/jeter1.bcf TMP/jeter2.vcf.gz
		countIt "gatkHardFilters" TMP/jeter1.vcf.gz TMP/jeter2.vcf.gz

		rm TMP/jeter2.vcf.gz TMP/jeter2.vcf.gz.tbi TMP/jeter1.vcf.gz TMP/jeter1.vcf.gz.tbi
		
	fi


	# min max alleles  ############################################################################
	bcftools view  -m2 -M ${max_alleles_count} -O u -o TMP/jeter2.bcf  TMP/jeter1.bcf
	countIt "too_many_alts" TMP/jeter1.bcf TMP/jeter2.bcf
	mv TMP/jeter2.bcf TMP/jeter1.bcf



	# remove in blaclisted regions ############################################################################
	bcftools view  --targets-overlap 2 --targets-file ^TMP/jeter.blacklisted.bed -O u -o TMP/jeter2.bcf TMP/jeter1.bcf
	countIt "blackListedRegions" TMP/jeter1.bcf TMP/jeter2.bcf
	mv TMP/jeter2.bcf TMP/jeter1.bcf
	

	# remove all annotations ################################################################################
	bcftools annotate --force \\
		-x '^INFO/AC,INFO/AN,INFO/ReadPosRankSum,INFO/MQRankSum,INFO/MQ,INFO/FS,INFO/QD,INFO/SOR,FILTER' \\
		-O u -o TMP/jeter2.bcf TMP/jeter1.bcf
	countIt "annotateX" TMP/jeter1.bcf TMP/jeter2.bcf
	mv TMP/jeter2.bcf TMP/jeter1.bcf


	# could happen that some variant are discarded here: saw some gatk4 variants where *NO* genotype was called. #############################
	if test -s TMP/cases.txt || test -s TMP/controls.txt ; then
		cat TMP/cases.txt TMP/controls.txt | sort -T TMP | uniq > TMP/jeter.samples
		bcftools view --trim-alt-alleles --samples-file 'TMP/jeter.samples' -O u -o TMP/jeter2.bcf TMP/jeter1.bcf
		rm TMP/jeter.samples
		countIt "samples" TMP/jeter1.bcf TMP/jeter2.bcf
		mv TMP/jeter2.bcf TMP/jeter1.bcf
	fi


	# force recalculation of AF/AC/AN #########################################################
	bcftools  +fill-tags -O u -o TMP/jeter2.bcf TMP/jeter1.bcf -- -t AN,AC,AF
	countIt "filltags" TMP/jeter1.bcf TMP/jeter2.bcf
	mv TMP/jeter2.bcf TMP/jeter1.bcf

	## too many no-call ("genotyping rate" >= 95%) #################################################
	bcftools view -i 'CHROM=="Y" || CHROM=="chrY" ||  F_MISSING < ${f_missing}' -O b -o TMP/jeter2.bcf TMP/jeter1.bcf
	countIt "F_MISSING_${f_missing}" TMP/jeter1.bcf TMP/jeter2.bcf
	mv TMP/jeter2.bcf TMP/jeter1.bcf

	## update VCF ID. Useful for plink stuff #########################################################
	if ${params.wgselect.with_setid == true } ; then
		bcftools annotate --set-id +'%VKX'  -O b -o TMP/jeter2.bcf TMP/jeter1.bcf
		mv TMP/jeter2.bcf TMP/jeter1.bcf
	fi

	
	## sex et homvar (1 homvar and 0 het) ################################################################
	if ${params.wgselect.with_homvar == true } ; then
		bcftools view -O v -o TMP/jeter2.vcf TMP/jeter1.bcf
		java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfpar TMP/jeter2.vcf > TMP/jeter1.vcf
		bcftools view -e 'INFO/SEX=0 &&  COUNT(GT="RA")==0 && COUNT(GT="AA")>0' -O u -o TMP/jeter2.bcf TMP/jeter1.vcf
		countIt "homvar" TMP/jeter1.bcf TMP/jeter2.bcf
		mv TMP/jeter2.bcf TMP/jeter1.bcf
		rm TMP/jeter1.vcf TMP/jeter2.vcf

	fi

	# split multiallelic #################################################################################""
	if [ "${max_alleles_count}" != "2" ] ; then

		bcftools norm -f "${reference}" --multiallelics -both  -O u -o TMP/jeter2.bcf TMP/jeter1.bcf
		countIt "norm" TMP/jeter1.bcf TMP/jeter2.bcf
		mv TMP/jeter2.bcf TMP/jeter1.bcf
	fi

	# not in pedigree ###########################################################################
	bcftools view -i 'AC[*]>0' -O v -o TMP/jeter2.vcf TMP/jeter1.bcf 

	countIt "AC_GT_0" TMP/jeter1.bcf TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf
	rm TMP/jeter1.bcf


	# ignore spanning deletions #################################################################
	awk -F '\t' '(\$0 ~ /^#/ || \$5!="*")'  TMP/jeter1.vcf > TMP/jeter2.vcf

	countIt "spandel" TMP/jeter1.vcf TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf


	## polyx ###################################################################################


	if ${ (polyx as int) > 1 } ; then
		echo "##show me vcf" 1>&2
		cat TMP/jeter1.vcf | head 1>&2

		java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar vcfpolyx -R "${reference}" --tag POLYX -n "${polyx}" TMP/jeter1.vcf |\
		bcftools view -e 'FILTER~"POLYX_ge_${polyx}"' > TMP/jeter2.vcf
		countIt "polyx${polyx}" TMP/jeter1.vcf TMP/jeter2.vcf

		mv TMP/jeter2.vcf TMP/jeter1.vcf
	fi
	
	## CADD ######################################################################################


	if ${!cadd_tabix.isEmpty() && (cadd_phred as double) > 0}  ; then
        	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar vcfcadd \
			--tabix "${cadd_tabix}" TMP/jeter1.vcf > TMP/jeter2.vcf
	      	mv TMP/jeter2.vcf TMP/jeter1.vcf
		
		bcftools view -e 'INFO/CADD_PHRED < ${cadd_phred as double}' -O v -o TMP/jeter2.vcf TMP/jeter1.vcf

		countIt "CADD" TMP/jeter1.vcf TMP/jeter2.vcf
	      	mv TMP/jeter2.vcf TMP/jeter1.vcf

	fi


	# test MAF
	
	if ${maxmaf>=0} && test -s TMP/pedigree.ped ; then
		java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfburdenmaf \
			--pedigree TMP/pedigree.ped --prefix "" --min-maf 0  --max-maf "${maxmaf}"  TMP/jeter1.vcf   > TMP/jeter2.vcf
		countIt "MAF" TMP/jeter1.vcf TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf

	fi

	# fisher per variant
	if ${fisherh >= 0.0} && test -s TMP/pedigree.ped ; then
		java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfburdenfisherh --filter '' \\
			--pedigree TMP/pedigree.ped \\
			--min-fisher "${fisherh}"  TMP/jeter1.vcf   > TMP/jeter2.vcf
		countIt "fisherH" TMP/jeter1.vcf TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf
		bcftools annotate -x 'FILTER/CTRL_CASE_RATIO' TMP/jeter1.vcf > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf

	fi

	# low or high DP

	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffilterjdk --nocode  -e 'final double dp= variant.getGenotypes().stream().filter(G->G.isCalled() && G.hasDP()).mapToInt(G->G.getDP()).average().orElse(${minDP}); return dp>=${minDP} && dp<=${maxDP};'  TMP/jeter1.vcf > TMP/jeter2.vcf
	countIt "lowDP" TMP/jeter1.vcf TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf


	# all genotypes with ALT must have GQ >= 'x'

	if [ "${lowGQ}" -gt 0 ] ; then
		## low GQ all genotypes carrying a ALT must be GQ > 'x'
		java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffilterjdk --nocode  -e 'return variant.getGenotypes().stream().filter(g->g.isCalled() && !g.isHomRef() && g.hasGQ()).allMatch(g->g.getGQ()>=${lowGQ});' TMP/jeter1.vcf > TMP/jeter2.vcf
		countIt "lowGQ${lowGQ}" TMP/jeter1.vcf TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf

	fi

	## singleton
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffilterjdk \\
		--nocode  \\
		-e 'Genotype singleton=null; for(final Genotype g: variant.getGenotypes()) {if(g.isCalled() && !g.isHomRef()) { if(singleton!=null) return true;singleton=g;}} if(singleton!=null && singleton.isFiltered()) return false; if(singleton!=null && singleton.isHet() && singleton.hasGQ() && singleton.getGQ()<${minGQsingleton}) return false; if(singleton !=null && singleton.hasAD() && singleton.isHet() && singleton.getAD().length==2) {int array[]=singleton.getAD();double r= array[1]/(double)(array[0]+array[1]);if(r< ${minRatioSingleton} || r>(1.0 - ${minRatioSingleton})) return false;} return true; ' \\
		TMP/jeter1.vcf > TMP/jeter2.vcf
	countIt "singleton" TMP/jeter1.vcf TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf

	if ${!mapability.isEmpty()} ; then

		# DukeMapability
		java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar vcfbigwig -tag mapability \
			-B "${mapability}" TMP/jeter1.vcf > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf

		java  -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffilterjdk --nocode   \
				-e 'final String tag= "mapability"; return !variant.hasAttribute(tag) ||  (variant.getAttributeAsDouble(tag,0.0)==1.0);' TMP/jeter1.vcf > TMP/jeter2.vcf
		countIt "mapability" TMP/jeter1.vcf TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf

	fi

	### GNOMAD GENOME #####################################################################################


	if ${params.wgselect.with_gnomad==true && !params.gnomad.isEmpty()} ; then

		test -f "${params.gnomad}"

        	# gnomad genome
        	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfgnomad \\
			--bufferSize 1000 \\
			-F '${gnomadPop}' \\
                	-g "${params.gnomad}" \\
			--max-af "${gnomadAF}" TMP/jeter1.vcf   > TMP/jeter2.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf

		if ${params.with_gnomad_filtered} ; then

			bcftools view --header-only TMP/jeter1.vcf |\\
				grep "^##FILTER" |\\
				tr "<," "\\n" |\\
				grep '^ID=GNOMAD_' |\\
				awk -F '=' '{printf("%s FILTER ~ \\"%s\\" ",(NR==1?"":" ||"),\$2);}' | tr '"' "'"  > TMP/gnomad.filters

			bcftools view -e "`cat TMP/gnomad.filters`"  -O v -o TMP/jeter2.vcf TMP/jeter1.vcf
			countIt "gnomadgenome.${gnomadPop}" TMP/jeter1.vcf TMP/jeter2.vcf

			mv TMP/jeter2.vcf TMP/jeter1.vcf
		fi	

	fi


 	## FUNCTIONNAL ANNOTATION ##############################################################################


	    if ${params.wgselect.with_vep==true} ; then

		vep_toto --output_file TMP/jeter2.vcf < TMP/jeter1.vcf
 		mv TMP/jeter2.vcf TMP/jeter1.vcf
	    fi


	    if  ${params.wgselect.with_snpeff==true} ; then 

	   	 # snpeff
		 snpEff -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  eff -dataDir "\${PWD}/${snpEffDir}" \\
			-interval "TMP/contig.bed" -nodownload -noNextProt -noMotif -noInteraction -noLog -noStats -chr chr -i vcf -o vcf "${params.snpeff_db}" TMP/jeter1.vcf > TMP/jeter2.vcf
	         mv TMP/jeter2.vcf TMP/jeter1.vcf

	    fi

	    if ${!soacn.isEmpty()} ; then
	    
	    	# filter prediction
		java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffilterso \
			${inverse_so?"--invert":""} \
			--remove-attribute  --rmnoatt \
			--acn "${soacn}" \
		   	TMP/jeter1.vcf  > TMP/jeter2.vcf

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

		countIt "exclude.prediction" TMP/jeter1.vcf TMP/jeter2.vcf
	        mv TMP/jeter2.vcf TMP/jeter1.vcf

	    fi


    # CONTRAST #############################################################################################

	if ${(params.wgselect.with_contrast as boolean)}  && test -s TMP/cases.txt && test -s TMP/controls.txt ; then
		bcftools +contrast \
			-0 "TMP/controls.txt" \
			-1 "TMP/cases.txt" \
			-a PASSOC,FASSOC,NASSOC,NOVELAL,NOVELGT -O v -o TMP/jeter2.vcf TMP/jeter1.vcf
		mv TMP/jeter2.vcf TMP/jeter1.vcf


	fi	


MD5=`cat TMP/jeter1.vcf | sha1sum | cut -d ' ' -f1`

bcftools view  -O u TMP/jeter1.vcf |\
bcftools sort -T TMP --max-mem "${task.memory.giga}G" -O b -o "\${MD5}.contig.bcf" 
bcftools index  "\${MD5}.contig.bcf"



mv TMP/contig.bed "\${MD5}.contig.bed"



countIt "final" TMP/jeter1.vcf contig.bcf
gzip --best TMP/variant_list.txt

mv TMP/variant_list.txt.gz "\${MD5}.variant_list.txt.gz"
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
	label "process_single"
	input:
        	path("DIR/*")
        output:
                path("wgselect.count.tsv"),emit:output
        script:
        	prefix = params.prefix?:""
        """
        hostname 1>&2
	set -o pipefail

        echo "#FILTER\tIN\tOUT\tDIFF" > "wgselect.count.tsv"
        find DIR -type l -name "*.gz" -exec gunzip -c '{}' ';'  |\\
		cut -f 1-4 |\\
		sort -T . -t '\t' -k1,1 |\\
        	datamash  -g 1  sum 2 sum 3 sum 4 >> wgselect.count.tsv
        	
        """
        }


