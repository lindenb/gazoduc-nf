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

include { moduleLoad; isHg19; isHg38} from '../../modules/utils/functions.nf'
include { BED_CLUSTER_01 } from '../../modules/jvarkit/jvarkit.bedcluster.01.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_CONCAT_01} from '../bcftools/bcftools.concat.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
//include {BED_WITH_VARIANTS_01} from '../../modules/bcftools/bed.with.variants.01.nf'

workflow SPLICEAI_01 {
	take:
		meta
		reference /* path to reference */
		vcf /* one indexed vcf or file with suffix .list with the path to the vcfs */
		pedigree /* jvarkit ped or NO_FILE */
		bed /* restrict to that bed or NO_FILE */
	main:
		version_ch = Channel.empty()
		exons_ch = DOWNLOAD_EXONS(meta,reference, bed)
		version_ch = version_ch.mix(exons_ch.version)

		//bed2_ch = BED_WITH_VARIANTS_01(meta, vcf.combine(exons_ch.bed))
		//version_ch = version_ch.mix(bed2_ch.version)

		clusters_ch = BED_CLUSTER_01(meta.plus("bed_cluster_method":"--jobs 1000"), reference, exons_ch.bed)
		version_ch = version_ch.mix(clusters_ch.version)
		
		each_bed = clusters_ch.output.splitText().map{T->file(T.trim())}

		spliceai_ch = APPLY_SPLICEAI(meta, reference, pedigree, vcf.combine(each_bed))
		version_ch = version_ch.mix(spliceai_ch.version)

		file_list_ch = COLLECT_TO_FILE_01([:],spliceai_ch.vcf.collect())
		version_ch = version_ch.mix(file_list_ch.version)

		concat_ch = BCFTOOLS_CONCAT_01([:],file_list_ch.output)
		version_ch = version_ch.mix(concat_ch.version)

                version_ch = MERGE_VERSION(meta, "SpliceAI", "SpliceAI", version_ch.collect())
	emit:
		version = version_ch
		vcf = concat_ch.vcf
		
	}

process DOWNLOAD_EXONS {
	afterScript "rm -rf TMP"
	input:
		val(meta)
		val(reference)
		path(bed)
	output:
		path("exons.bed.gz"),emit:bed
		path("version.xml"),emit:version
	script:
		def distance = meta.splice_distance?:50
		def url = isHg19(reference)?"https://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeBasicV24lift37.txt.gz":""
	if(isHg19(reference))
	"""
	hostname 1>&2
	${moduleLoad("jvarkit bedtools")}
	set -o pipefail

	mkdir TMP

	if ${!bed.name.equals("NO_FILE")} ; then
		java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP "${bed}" |\
			sort -T TMP  -t '\t' -k1,1 -k2,2n |\
                	bedtools merge > TMP/jeter2.bed
	fi


	# url dans le help de spliceai
	wget -O - "${url}" |\
		gunzip -c |\
		awk -F '\t' '{N=int(\$9);X=${distance};split(\$10,S,/[,]/);split(\$11,E,/[,]/);for(i=1;i<=N;i++) {V=int(S[i]);if(i>1) printf("%s\t%d\t%d\\n",\$3,(V<X?0:V-X),V+X);V=int(E[i]);if(i<N) printf("%s\t%d\t%d\\n",\$3,(V<X?0:V-X),V+X);} }' |\
		java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP |\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\
		bedtools merge |\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\
		${bed.name.equals("NO_FILE")?"":"bedtools intersect -u -a - -b TMP/jeter2.bed |"} \
		gzip --best > TMP/exons.bed.gz

	mv TMP/exons.bed.gz ./


	##################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">extract splice interval as bed file</entry>
		<entry key="url"><url>${url}</url></entry>
	        <entry key="splice_distance">${distance}</entry>
	</properties>
	EOF
	"""
	}

process APPLY_SPLICEAI {
	tag "${bed.name}"
	conda "${params.conda}/SPLICEAI"
	afterScript "rm -rf TMP"
	cpus 2
	input:
		val(meta)
		val(reference)
		path(pedigree)
		tuple path(vcf),path(bed)
	output:
		path("contig.bcf"),emit:vcf
		path("contig.bcf.csi"),emit:index
		path("version.xml"),emit:version
	script:
		def db = isHg19(reference)?"grch37":"UNDEFINED"
		def distance = meta.splice_distance?:50
	"""
	hostname 1>&2
	module load bcftools/0.0.0
	set -o pipefail

	export OMP_NUM_THREADS=${task.cpus}

	mkdir TMP
	export TMPDIR=\${PWD}/TMP
		
	if ${vcf.name.endsWith(".list")} ; then	
		bcftools concat --remove-duplicates  --allow-overlaps --file-list "${vcf.toRealPath()}" --regions-file "${bed}" -O b -o TMP/jeter.bcf
	else
		bcftools view --regions-file "${bed}" -O b -o TMP/jeter.bcf "${vcf.toRealPath()}"
	fi

	if ${!pedigree.name.equals("NO_FILE")} ; then
 
		grep -v '^#' "${pedigree}" | cut -f 2 | sort -T TMP | uniq > TMP/a
		bcftools query -l TMP/jeter.bcf  | sort -T TMP | uniq > TMP/b
		comm -12 TMP/a TMP/b > TMP/jeter.samples
		test -s TMP/jeter.samples

		awk -F '\t' '(\$6=="affected" || \$6=="case")'  "${pedigree}" | cut -f 2 | sort -T TMP | uniq > TMP/cases1.txt
		awk -F '\t' '(\$6=="unaffected" || \$6=="control")'  "${pedigree}" | cut -f 2  | sort -T TMP | uniq > TMP/controls1.txt
		comm -12 TMP/b TMP/cases1.txt > TMP/cases.txt
		comm -12 TMP/b TMP/controls1.txt > TMP/controls.txt

		bcftools view --samples-file TMP/jeter.samples -O b -o TMP/jeter2.bcf TMP/jeter.bcf
		mv TMP/jeter2.bcf TMP/jeter.bcf
	fi

	bcftools annotate -O u -x '^INFO/AC,^INFO/AF,FILTER' TMP/jeter.bcf  |\
	        bcftools norm -f "${reference}" --multiallelics -both -O u |\
		bcftools view --min-ac 1 -e 'ALT="*" ' -O b -o TMP/jeter2.bcf
	mv TMP/jeter2.bcf TMP/jeter.bcf

	bcftools view TMP/jeter.bcf |\
	        spliceai -R "${reference}"  -A "${db}" -D ${distance} |\
        	bcftools view  -i 'INFO/SpliceAI!=""' -O b -o TMP/jeter2.bcf
	mv TMP/jeter2.bcf TMP/jeter.bcf


	if test -s TMP/cases.txt && test -s TMP/controls.txt ; then
		
		bcftools +contrast \
		     -0 "TMP/controls.txt" \
                     -1 "TMP/cases.txt" \
		     -a PASSOC,FASSOC,NASSOC,NOVELAL,NOVELGT \
		     -O b -o TMP/jeter2.bcf TMP/jeter.bcf

		mv TMP/jeter2.bcf TMP/jeter.bcf
	fi


	
	bcftools index TMP/jeter.bcf

	mv TMP/jeter.bcf ./contig.bcf 
	mv TMP/jeter.bcf.csi ./contig.bcf.csi

	##################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">apply spliceai</entry>
		<entry key="db">${db}</entry>
		<entry key="pedigree">${pedigree}</entry>
		<entry key="bed">${bed}</entry>
	        <entry key="splice_distance">${distance}</entry>
		<entry key="bcftools.version">\$( bcftools --version-only)</entry>
		<entry key="spliceai.version">\$(spliceai --help | grep Version)</entry>
	</properties>
	EOF
	"""
	}
