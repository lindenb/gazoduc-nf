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



include {getVersionCmd;runOnComplete;moduleLoad;parseBoolean;dumpParams;paramsToString} from '../../../modules/utils/functions.nf'
include {VALIDATE_CASE_CONTROL_PED_01} from '../../../modules/pedigree/validate.case.ctrl.pedigree.01.nf'
include {VCF_INTER_CASES_CONTROLS_01} from '../../../subworkflows/bcftools/vcf.inter.cases.controls.01.nf'
include {BED_CLUSTER_01} from '../../../modules/jvarkit/jvarkit.bedcluster.01.nf'
include {SQRT_FILE} from '../../../modules/utils/sqrt.nf'
include {COLLECT_TO_FILE_01} from '../../../modules/utils/collect2file.01.nf'
include {WGSELECT_01} from '../../../subworkflows/wgselect/wgselect.01.nf'
include {PEDIGREE_FOR_RVTESTS} from '../../../modules/rvtests/rvtests.cases.controls.ped.01.nf'
include {RVTESTS_POST_PROCESS} from '../../../subworkflows/rvtests/rvtests.post.process.01.nf'
include {CONCAT_FILES_01} from '../../../modules/utils/concat.files.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {MULTIQC_01} from '../../../modules/multiqc/multiqc.01.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.02.nf'
include {PIHAT_CASES_CONTROLS_01} from '../../../subworkflows/pihat/pihat.cases.controls.01.nf'
include {VCF_TO_BED} from '../../../modules/bcftools/vcf2bed.01.nf'
include {JVARKIT_GATK_HARD_FILTERING_01} from '../../../subworkflows/jvarkit/jvarkit.gatk_hard_filtering.01.nf'
include {PARAMS_MULTIQC} from '../../../modules/utils/params.multiqc.nf'

if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}




workflow {
		BURDEN_CODING(params.genomeId, Channel.fromPath(params.vcf), file(params.pedigree), file(params.bed))
		}

workflow BURDEN_CODING {
	take:
		genomeId
		vcf
		pedigree
		bed
	main:
		to_zip = Channel.empty()
		version_ch = Channel.empty()


		
		ped_ch = VALIDATE_CASE_CONTROL_PED_01(pedigree)
		version_ch = version_ch.mix(ped_ch.version)

		vcf_inter_ch = VCF_INTER_CASES_CONTROLS_01(
				vcf,
				ped_ch.cases_list, ped_ch.controls_list
				)
		version_ch = version_ch.mix(vcf_inter_ch.version)

		if(parseBoolean(params.burden.with_pihat)) {
			pihat = PIHAT_CASES_CONTROLS_01(genomeId,file(vcf),ped_ch.cases_list,ped_ch.controls_list)
			version_ch = version_ch.mix(pihat.version)
			to_zip = to_zip.mix(pihat.pihat_png)
			to_zip = to_zip.mix(pihat.pihat_sample2avg_png)
			to_zip = to_zip.mix(pihat.removed_samples)
			to_zip = to_zip.mix(pihat.plink_genome)

			rebuild_ch = REBUILD_PEDIGREE(pedigree,pihat.cases,pihat.controls)
			version_ch = version_ch.mix(rebuild_ch.version)

			new_ped_ch = rebuild_ch.pedigree
			}
		else
			{
			rebuild_ch = REBUILD_PEDIGREE(pedigree,ped_ch.cases_list,ped_ch.controls_list)
			version_ch = version_ch.mix(rebuild_ch.version)

			new_ped_ch = rebuild_ch.pedigree
			}

		vcf2bed_ch = VCF_TO_BED([with_header: false],vcf)
		version_ch = version_ch.mix(vcf2bed_ch.version)

		/* compute GATK Hard filters if needed */

		if((params.wgselect.gatk_hardfiltering_percentile as int)  > 0 ) {
			gatk_hard_filters_ch = JVARKIT_GATK_HARD_FILTERING_01(
				[percentile: params.wgselect.gatk_hardfiltering_percentile],
				vcf2bed_ch.bed.splitCsv(sep:'\t', header:false).map{T->[
					interval: T[0]+":"+((T[1] as int)+1)+"-"+T[2],
					vcf: T[3]
					]}
				)
			version_ch = version_ch.mix(gatk_hard_filters_ch.version)
			hard_filters_ch = gatk_hard_filters_ch.output
			}
		else  
			{
			hard_filters_ch = Channel.fromPath(file("NO_FILE"))
			}
		


		vcfinterbed_ch = INTERSECT_VCF_WITH_USER_BED(vcf2bed_ch.bed,bed)
		version_ch = version_ch.mix(vcfinterbed_ch.version)
		
		genes_ch = EXTRACT_GENES(genomeId, vcf2bed_ch.bed)
		version_ch = version_ch.mix(genes_ch.version)

		exons_ch =  EXTRACT_EXONS(genomeId, genes_ch.bed)
		version_ch = version_ch.mix(exons_ch.version)

		cluster_ch = BED_CLUSTER_01(
                        [bed_cluster_method : params.bed_cluster_method],
                        genomeId,
                        exons_ch.bed
                        )
                version_ch = version_ch.mix(cluster_ch.version)


		wgselect_ch = WGSELECT_01(
				[:],
				genomeId, 
				cluster_ch.output.splitText().map{it.trim()}.
				combine(vcf).
				combine(new_ped_ch).
				combine(hard_filters_ch).
				map{T->[
					"bed": file(T[0]),
					"vcf": file(T[1]),
					"pedigree" : T[2],
					"hard_filters": T[3]
				]})


		version_ch = version_ch.mix(wgselect_ch.version)
		to_zip = to_zip.mix(wgselect_ch.variants_list)


		rvtests_ped_ch = PEDIGREE_FOR_RVTESTS([:], new_ped_ch)
		version_ch = version_ch.mix(rvtests_ped_ch.version)



		per_gene_ch = RVTEST_PER_GENES(
			[:],
			genomeId,
			rvtests_ped_ch.pedigree,
			wgselect_ch.bed,
			genes_ch.bed.splitCsv(header:true,sep:'\t').collate(params.collate_size)
			)
		version_ch = version_ch.mix(per_gene_ch.version)


		concat_ch = CONCAT_FILES_01([suffix:".list",concat_n_files:50,downstream_cmd:""], per_gene_ch.output.collect())
		version_ch = version_ch.mix(concat_ch.version)


		digest_ch = RVTESTS_POST_PROCESS([:],genomeId, "${params.prefix?:""}BurdenCoding" ,concat_ch.output)
                version_ch = version_ch.mix(digest_ch.version)
		to_zip = to_zip.mix(digest_ch.zip)

		params4multiqc_ch = PARAMS_MULTIQC([:])
		


		multiqc_ch = MULTIQC_01(["title":"${params.prefix}BurdenCoding","comment":"VCF ${params.vcf}"], digest_ch.to_multiqc.concat(params4multiqc_ch.output).collect());
                version_ch = version_ch.mix(multiqc_ch.version)

		version_ch = MERGE_VERSION("burden coding",version_ch.collect())

		html = VERSION_TO_HTML(version_ch.version)
		to_zip = to_zip.mix(html.html)


		
		PUBLISH(to_zip.collect())

	}


process INTERSECT_VCF_WITH_USER_BED {
tag "${vcf2bed.name} ${userbed.name}"
input:
	path(vcf2bed)
	path(userbed)
output:
	path("bedintervcf.bed"),emit:bed
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("htslib jvarkit bedtools")}
set -o pipefail
mkdir -p TMP

if ${userbed.name.equals("NO_FILE")} ; then

	cut -f1,2,3 '${vcf2bed}'  | LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n  | bedtools merge > bedintervcf.bed

else

	cut -f1,2,3 '${vcf2bed}'  | LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n  | bedtools merge > TMP/a.bed
	cut -f1,2,3 '${userbed}'  | LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n  | bedtools merge > TMP/b.bed

	bedtools intersect -a TMP/a.bed -b TMP/b.bed | LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n > bedintervcf.bed

fi

test -s bedintervcf.bed

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
</properties>
EOF
"""
}         


process EXTRACT_GENES {
memory "2g"
afterScript "rm -rf TMP"
input:
	val(genomeId)
	path(bed) //bed is already filtered for intersect with vcf
output:
	path("genes.bed"),emit:bed
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
	def gtf = genome.gtf
	def slop= params.genes_slop
"""
hostname 1>&2
${moduleLoad("htslib jvarkit bedtools")}
set -o pipefail
mkdir -p TMP


cat << EOF | paste -s -d '\t' > TMP/tmp1.bed
chrom
start
end
feature
gene_biotype
gene_name
gene_id
EOF

tabix --regions "${bed}" "${gtf}" |\
	java -jar \${JVARKIT_DIST}/gtf2bed.jar -R "${reference}" --columns "gtf.feature,gene_biotype,gene_name,gene_id" |\
	awk -F '\t' '(\$4=="gene" && \$5=="protein_coding" && \$5!="." && \$6!=".")' |\
	bedtools slop -i -  -g "${reference}.fai" -b ${slop} |\
	LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n |\
	bedtools intersect -u -wa -a - -b ${bed} |\
	LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n > TMP/tmp2.bed

test -s TMP/tmp2.bed

cat TMP/tmp1.bed TMP/tmp2.bed > genes.bed

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">extract protein_coding genes from gtf</entry>
	<entry key="gtf">${gtf}</entry>
	<entry key="versions">${getVersionCmd("bedtools awk jvarkit/gtf2bed tabix")}</entry>
	<entry key="bed">${bed}</entry>
</properties>
EOF
"""
}


process EXTRACT_EXONS {
memory "2g"
afterScript "rm -rf TMP"
input:
	val(genomeId)
	path(genesbed)
output:
	path("exons.bed"),emit:bed
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
	def gtf = genome.gtf
	def slop= params.genes_slop
"""
hostname 1>&2
${moduleLoad("htslib jvarkit bedtools")}
set -o pipefail
mkdir -p TMP

# remove header
awk -F '\t' '(\$1!="chrom")' '${genesbed}' > TMP/genes.bed

tabix --regions TMP/genes.bed  "${gtf}" |\
	java -jar \${JVARKIT_DIST}/gtf2bed.jar -R "${reference}" --columns "gtf.feature"|\
	awk -F '\t' '(\$4=="exon")' |\
	cut -f1,2,3 |\
	bedtools slop -i - -g "${reference}.fai" -b ${slop} |\
	LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n |\
	bedtools merge > TMP/tmp1.bed

mv TMP/tmp1.bed exons.bed

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">extract protein_coding genes from gtf</entry>
	<entry key="gtf">${gtf}</entry>
	<entry key="versions">${getVersionCmd("bedtools awk jvarkit/gtf2bed tabix")}</entry>
	<entry key="bed">${genesbed}</entry>
</properties>
EOF
"""
}



process RVTEST_PER_GENES {
tag "${L.collect{T->T.gene_name}.join(",")}"
afterScript "rm -rf TMP"
memory {task.attempt==1?'5G':'20G'}
errorStrategy 'retry'
maxRetries 3
input:
	val(meta)
	val(genomeId)
	path(pedigree)
	path(vcfbed)
	val(L)
errorStrategy "retry"

output:
	path("assoc.list"),emit:output
	path("version.xml"),emit:version
script:
	if(!params.burden.containsKey("rvtest_arguments")) throw new IllegalArgumentException("params.burden.rvtest_arguments missing");
	def rvtest_arguments = params.burden.rvtest_arguments
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
"""
hostname 1>&2
${moduleLoad("rvtests bcftools jvarkit bedtools")}

mkdir -p TMP/ASSOC ASSOC

i=1


## https://github.com/zhanxw/rvtests/issues/80
cut -f 1 "${reference}.fai" > TMP/chroms.A.txt
sed 's/^chr//' TMP/chroms.A.txt > TMP/chroms.B.txt
paste TMP/chroms.A.txt TMP/chroms.B.txt > TMP/chroms.C.txt


# sort VCF bed
sort -t '\t' -T TMP -k1,1 -k2,2n  '${vcfbed}' > TMP/vcfs.bed


cat << EOF > TMP/input.tsv
${L.collect(T->T.chrom+"\t"+T.start+"\t"+T.end+"\t"+T.gene_name+"\t"+T.gene_id).join("\n")}
EOF

cat TMP/input.tsv | while IFS=\$'\\t' read CHROM START END GENE_NAME GENE_ID
do

# create a bed for this gene
echo "\${CHROM}\t\${START}\t\${END}" > TMP/gene.bed


# find VCF overlapping this gene
bedtools intersect -wa -u -a TMP/vcfs.bed -b TMP/gene.bed | cut -f4 | sort | uniq > TMP/vcfs.list


# vcf is empty (wgselect removed al the genes) . peek a random one
if test ! -s TMP/vcfs.list ; then
	head -n 1 TMP/vcfs.bed | cut -f 4 > TMP/vcfs.list
fi
test -s TMP/vcfs.list

# file containing the name of the gene
echo "\${GENE_ID}" > TMP/gene.id



bcftools concat --regions-file TMP/gene.bed --file-list TMP/vcfs.list -O u --allow-overlaps --remove-duplicates |\
	bcftools annotate -O v --rename-chrs TMP/chroms.C.txt |\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffiltergenes -g TMP/gene.id > TMP/jeter.vcf

	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfgenesplitter \
		-E 'ANN/GeneId ANN/FeatureId VEP/GeneId VEP/Ensp VEP/Feature' \
		--manifest TMP/jeter.mf \
		-o \${PWD}/TMP \
		TMP/jeter.vcf

	rm TMP/jeter.vcf

	awk -F '\t' '/^[^#]/ {S=\$4;gsub(/[\\/]+/,"_",S);G=\$5;gsub(/[\\/_\\-]+/,"_",G);K=\$6;gsub(/[\\/]+/,"_",K);printf("%s_%s_%s\t%s:%d-%d\t%s\\n",K,G,S,\$1,\$2,\$3,\$7);}' TMP/jeter.mf |\
	while read K RANGE V
	do
		bcftools index -t "\${V}"

		# build setFile
		echo "\$K\t\${RANGE}" > TMP/variants.setfile

		rvtest  --noweb \
        		--inVcf "\${V}" \
			--setFile TMP/variants.setfile \
	        	--pheno "${pedigree}" \
		        --out "TMP/ASSOC/part.\${GENE_ID}.0\${i}" \
			${rvtest_arguments} 2> TMP/last.rvtest.log

		rm "\${V}" "\${V}.tbi"

		i=\$((i+1))
	done

done

# collect each assoc by type and merge
find \${PWD}/TMP/ASSOC -type f -name "part.*assoc" | awk -F '/' '{N=split(\$NF,a,/\\./); print a[N-1];}' | sort -T TMP | uniq | while read SUFFIX
do

	find \${PWD}/TMP/ASSOC -type f -name "part.*\${SUFFIX}.assoc" > TMP/assoc.list
	test -s TMP/assoc.list

	xargs -a TMP/assoc.list cat | sed 's/^Range/#Range/' | LC_ALL=C sort -T TMP | uniq | sed 's/^#Range/Range/' > "ASSOC/concat.\${SUFFIX}.assoc"
done


find \${PWD}/ASSOC -type f -name "concat.*.assoc" >  assoc.list


cat << EOF > version.xml
<properties id="${task.process}">
  <entry key="name">${task.process}</entry>
  <entry key="description">VCF is split by transcript / gene, and then rvtest is applied.</entry>
  <entry key="name"Pedigree">${pedigree}</entry>
  <entry key="rvtest.path">\$(which rvtest)</entry>
  <entry key="rvtest.version">\$(rvtest  --version 2> /dev/null  | head -n1)</entry>
</properties>
EOF
"""
stub:
"""
touch assoc.list
echo "<properties/>" > version.xml
"""
}


process REBUILD_PEDIGREE {
executor "local"
input:
	path(pedigree)
	path(cases)
	path(controls)
output:
	path("case_control.ped"),emit:pedigree
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
set -o pipefail

sort -T . -t '\t' -k2,2 "${pedigree}" > tmp1.ped
cat  "${cases}" "${controls}" | sort | uniq > tmp2.ped

join -t '\t' -1 2 -2 1 -o '1.1,1.2,1.3,1.4,1.5,1.6' tmp1.ped  tmp2.ped |\
	sort -T . -t '\t' -k1,1 -k2,2  > case_control.ped
test -s case_control.ped
rm tmp1.ped tmp2.ped



###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">make pedigree after samples have been removed.</entry>
	<entry key="pedigree">${pedigree.toRealPath()}</entry>
	<entry key="cases">${cases.toRealPath()}</entry>
	<entry key="controls">${controls.toRealPath()}</entry>
	<entry key="samples.removed">\$(comm -23 <(cut -f 2 "${pedigree}" | sort)   <(cut -f 2 "case_control.ped" | sort)  | paste -s -d ' ' )</entry>
</properties>
EOF
"""
}



process PUBLISH {
tag "N=${files.size()}"
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
        val(files)
output:
        path("${params.prefix?:""}archive.zip")
when:
        !params.getOrDefault("publishDir","").trim().isEmpty()
script:
        prefix = params.getOrDefault("prefix","")
"""

mkdir "${prefix}archive"

cat << EOF | while read F ; do ln -s "\${F}" "./${prefix}archive/" ; done
${files.join("\n")}
EOF

zip -r "${prefix}archive.zip" "${prefix}archive"
"""
}

runOnComplete(workflow);
