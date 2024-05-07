/*

Copyright (c) 2024 Pierre Lindenbaum

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

include {slurpJsonFile;moduleLoad} from '../../modules/utils/functions.nf'
include {hasFeature;isBlank;backDelete} from './annot.functions.nf'


def TAG="SPLICEAI"

 
workflow ANNOTATE_SPLICEAI {
	take:
		genomeId
		vcfs /** json vcf,vcf_index */
	main:
		if(hasFeature("spliceai") && !isBlank(params.genomes[genomeId],"gtf") && !isBlank(params.genomes[genomeId],"spliceai_annotation_type")) {
			annotate_ch = ANNOTATE(genomeId,vcfs)
			out1 = annotate_ch.output
			out2 = annotate_ch.count
			out3 = MAKE_DOC(genomeId).output

			SAVE_SCORES(genomeId,annotate_ch.scores.collect())
			}
		else
			{
			out1 = vcfs
			out2 = Channel.empty()
			out3 = Channel.empty()
			}
	emit:
		output = out1
		count = out2
		doc = out3
}


process MAKE_DOC {
executor "local"
input:
        val(genomeId)
output:
	path("${TAG}.html"),emit:output
script:
	def genome = params.genomes[genomeId]
"""
cat << __EOF__ > ${TAG}.html
<dl>
<dt>${TAG}</dt>
<dd>SpliceAI Prediction</dd>
</dl>
__EOF__
"""
}



workflow ANNOTATE {
	take:
		genomeId
		json
	main:
		ch1 = SPLIT(genomeId, json)

		/** https://nextflow.slack.com/archives/C02T98A23U7/p1709801086681969 */
		ch1a = ch1.vcfs.
			map{it[1] instanceof List?it:[it[0],Collections.singletonList(it[1])]}.
			flatMap(KV->KV[1].collect{VCF->[KV[0],VCF]})

		ch1b = ch1.tbis.
			map{it[1] instanceof List?it:[it[0],Collections.singletonList(it[1])]}.
			flatMap(KV->KV[1].collect{TBI->[KV[0],TBI]})

		/*
		ch1a.view{"#DEBUGA "+ it+"\n"}
		ch1b.view{"#DEBUGB "+ it+"\n"}
		ch1c = ch1a.join(ch1b,failOnMismatch:true).view{"#DEBUGC "+ it+"\n"}
		ch1c = Channel.empty()
		*/
		
		ch1c = ch1a.join(ch1b,failOnMismatch:true , failOnDuplicate:false /* ah ? */)
		
		ch2 = APPLY_SPLICEAI(genomeId,ch1c )
		
		ch3 = JOIN(ch1.off.
				mix(ch2.output).
				groupTuple().
				map{T->[T[0], T[1].sort(), T[2].sort()]} 
			)
	emit:
		
		output = ch3.output
		count  = ch3.count
		scores = ch3.scores
	}


process SPLIT {
tag "${json.name}"
afterScript "rm -rf TMP"
memory "3g"
conda "/CONDAS/users/lindenbaum-p/SPLICEAI"
input:
	val(genomeId)
	//tuple path(vcf),path(vcf_idx),path(bed)
	path(json)
output:
	tuple path(json),path("OUTPUT/CHUNCKS/*.vcf.gz"),emit:vcfs
	tuple path(json),path("OUTPUT/CHUNCKS/*.vcf.gz.tbi"),emit:tbis
	tuple path(json),path("OUTPUT/off.bcf"),path("OUTPUT/off.bcf.csi"),emit:off
script:
	def genome = params.genomes[genomeId]
	def gtf = genome.gtf
	def distance = params.annotations.spliceai_distance?:50
	def count = 1000
	def row = slurpJsonFile(json)
	def ac = params.annotations.spliceai_max_ac?:0
"""
hostname 1>&2
${moduleLoad("bcftools bedtools htslib jvarkit")}
mkdir -p TMP/CHUNCKS


tabix --regions "${row.bed}" "${gtf}" |\\
	awk -F '\t' '(\$3=="exon") {for(i=0;i<2;i++) {P=(i==0?int(\$4)-1:int(\$5));printf("%s\t%d\t%d\\n",\$1,P,P+1);}}' |\
	bedtools slop -b "${distance}" -g "${genome.fasta}.fai" |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\
	bedtools merge |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n  > TMP/junctions.bed
	
if ! test -s TMP/junctions.bed
then
	tail -n1 "${genome.fasta}.fai"  | awk -F '\t' '{printf("%s\t0\t1\\n",\$1);}' >  TMP/junctions.bed
fi

bgzip --force TMP/junctions.bed
tabix --force -p bed TMP/junctions.bed.gz

bcftools view -O b --targets-file ^TMP/junctions.bed.gz -o TMP/off1.bcf  "${row.vcf}"
bcftools index --force TMP/off1.bcf 

bcftools view -O b --targets-file TMP/junctions.bed.gz --min-ac '${(ac as int)+1}:nref'  -o TMP/off2.bcf "${row.vcf}"
bcftools index --force TMP/off2.bcf 

bcftools concat --remove-duplicates --allow-overlaps -O b -o  TMP/off.bcf TMP/off1.bcf TMP/off2.bcf
bcftools index --force TMP/off.bcf

bcftools view --targets-file TMP/junctions.bed.gz "${row.vcf}" --max-ac '${ac}:nref' |\\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar vcfsplitnvariants --variants-count ${count} --index -o TMP/CHUNCKS/split

# create one if no VCF above
ls TMP/CHUNCKS/split.*.vcf.gz > /dev/null || ( bcftools view --header-only -O z -o TMP/CHUNCKS/split.0.vcf.gz '${row.vcf}' && bcftools index -t TMP/CHUNCKS/split.0.vcf.gz )

find . -type f -name "*.bcf" -o -name "*.vcf.gz" |while read F; do echo -n "\${F}: " 1>&2 && bcftools query -N -f '.' "\${F}" | wc -c 1>&2 ; done

mkdir -p OUTPUT
mv TMP/off.bcf OUTPUT/
mv TMP/off.bcf.csi OUTPUT/
mv TMP/CHUNCKS OUTPUT/

"""
}

process APPLY_SPLICEAI {
tag "${vcf.name} ${json.name}"
afterScript "rm -rf TMP"
cpus 5
conda "/CONDAS/users/lindenbaum-p/SPLICEAI"
input:
	val(genomeId)
	tuple path(json),path(vcf),path(tbi)
output:
	tuple path(json),path("OUTPUT/${TAG}.bcf"),path("OUTPUT/${TAG}.bcf.csi"),emit:output
script:
	def genome = params.genomes[genomeId]
	def db =genome.spliceai_annotation_type
	def distance = params.annotations.spliceai_distance?:50
	def treshold = params.annotations.spliceai_highscore?:0.9
	def ac = params.annotations.spliceai_max_ac?:0
"""
hostname 1>&2
${moduleLoad("bcftools htslib")}
mkdir -p TMP OUTPUT

echo "\${PATH}"
which conda 1>&2 | cat

set -x


# previous scores where definied
if ${!file(params.annotations.spliceai_tabix_scores).name.equals("NO_FILE")}
then

	echo '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3.1 variant annotation. These include delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">' > TMP/tmp.header

	bcftools annotate -a ${params.annotations.spliceai_tabix_scores}' \\
			-h TMP/tmp.header \\
			-c "CHROM,POS,REF,ALT,SpliceAI" \\
			-O b -o TMP/tmp1.bcf '${vcf}'


	# splice AI was set exclude flag
	bcftools view -e 'SpliceAI==""' -O b -o TMP/done.bcf TMP/tmp1.bcf
	bcftools index TMP/done.bcf


	# splice AI was NOT set
	bcftools view -i 'SpliceAI==""' -O z -o TMP/remain.vcf.gz TMP/tmp1.bcf
	bcftools index -t TMP/remain.vcf.gz

	rm TMP/tmp1.bcf

else
	# empty vcf
	bcftools view --header-only -O b -o TMP/done.bcf '${vcf}'
	bcftools index TMP/done.bcf

	# copy original vcf
	bcftools view --threads ${tas.cpu} -o TMP/remain.vcf.gz -O z '${vcf}'
	bcftools index -t TMP/remain.vcf.gz
fi


export OMP_NUM_THREADS=${task.cpus}

spliceai -R "${genome.fasta}"  -A "${db}" -D ${distance} -I TMP/remain.vcf.gz  |\\
	bcftools view -O b -o TMP/in.bcf

# merge with done vcf
bcftools concat -a -O b -o TMP/tmp1.bcf TMP/in.bcf TMP/done.bcf
mv TMP/tmp1.bcf TMP/in.bcf

# search high spliceai scores
bcftools  query -f '%CHROM\t%POS\t%REF\t%ALT\t%SpliceAI\\n'  TMP/in.bcf |\\
	LC_ALL=C awk -F '\t' '{n1=split(\$5,a,/,/);for(i=1;i<=n1;i++) {split(a[i],b,/\\|/);if (b[3] >= ${treshold} || b[4] >= ${treshold} || b[5] >= ${treshold} || b[6]>= ${treshold}) {printf("%s\t%s\t%s\t%s\t1\\n",\$1,\$2,\$3,\$4); break;}}}' |\\
	LC_ALL=C sort -t '\t' -k1,1 -k2,2n -T TMP |\\
	bgzip > TMP/high.tsv.gz

gunzip -c TMP/high.tsv.gz | head 1>&2

tabix --force -s 1 -b 2 -e 2 TMP/high.tsv.gz

echo '##INFO=<ID=${TAG}_HIGH,Number=0,Type=Flag,Description="SpliceAI score>=${treshold} evaluated only when INFO/AC <= ${ac} .">' > TMP/high.header

# annotate high scores with FLAG
bcftools annotate -a TMP/high.tsv.gz -h TMP/high.header -c "CHROM,POS,REF,ALT,${TAG}_HIGH" -O b -o TMP/in2.bcf TMP/in.bcf
mv TMP/in2.bcf TMP/${TAG}.bcf
	
bcftools index --force TMP/${TAG}.bcf
mv TMP/${TAG}.* OUTPUT/
"""
}


process JOIN {
tag "${json.name} N=${vcfs}+${vcfs_idx}"
afterScript "rm -rf TMP"
input:
	tuple path(json),path("dir??/*"),path("dir??/*")
output:
	path("OUTPUT/${TAG}.json"),emit:output
	path("OUTPUT/${TAG}.count"),emit:count
	path("OUTPUT/${TAG}.scores.tsv"),emit:scores
script:
	def row = slurpJsonFile(json)
"""
hostname 1>&2
${moduleLoad("bcftools bedtools htslib")}
mkdir -p TMP OUTPUT

find ./dir* -name "*.vcf.gz" -o -name "*.bcf" > TMP/jeter.list

bcftools concat --allow-overlaps -O b -o  TMP/${TAG}.bcf --file-list TMP/jeter.list
bcftools index --force TMP/${TAG}.bcf

# save spliceAI score for another day
bcftools query -i 'INFO/SpliceAI!=""' -f '%CHROM\t%POS\t%REF\t%SpliceAI\\n' TMP/${TAG}.bcf |\\
	awk -F '\t' '{n=split(\$4,a,/[,]/);for(i=1;i<=n;i++)  {split(a[i],b,/[\|]/); printf("%s\t%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,b[1],a[i]);}}' |\\
	LC_ALL=C sort -T TMP -t \$'\t' -k1,1 -k2,2n -k3,3 -k4,4 > TMP/${TAG}.scores.tsv



cat << EOF > TMP/${TAG}.json
{
"vcf"   : "\${PWD}/OUTPUT/${TAG}.bcf",
"index" : "\${PWD}/OUTPUT/${TAG}.bcf.csi",
"bed"   : "${row.bed}"
}
EOF


bcftools query -N -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/${TAG}.count
mv -v TMP/${TAG}.* OUTPUT/
${backDelete(row)}
"""
}


process SAVE_SCORES {
tag "N=${L.size()}"
afterScript "rm -rf TMP"
input:
	val(genomeId)
	val(L)
output:
	tuple path("${params.prefix?:""}${genomeId}.spliceAI.scores.tsv.gz"),path("${params.prefix?:""}${genomeId}.spliceAI.scores.tsv.gz.tbi"),emit:output
script:
"""
hostname 1>&2
${moduleLoad("htslib")}

mkdir -p TMP
LC_ALL=C sort -T TMP -t \$'\t' -k1,1 -k2,2n -k3,3 -k4,4 ${L.join(" ")} --merge -o "TMP/${params.prefix?:""}${genomeId}.spliceAI.scores.tsv"


if ${!file(params.annotations.spliceai_tabix_scores).name.equals("NO_FILE")}
then
	gunzip -c "${params.annotations.spliceai_tabix_scores}" > TMP/prev.tsv
	LC_ALL=C sort -T TMP -t \$'\t' -k1,1 -k2,2n -k3,3 -k4,4 --merge  "TMP/${params.prefix?:""}${genomeId}.spliceAI.scores.tsv" TMP/prev.tsv -o "TMP/${params.prefix?:""}${genomeId}.spliceAI.scores.merged.tsv
	rm TMP/prev.tsv
else
	cp "TMP/${params.prefix?:""}${genomeId}.spliceAI.scores.tsv"  "TMP/${params.prefix?:""}${genomeId}.spliceAI.scores.merged.tsv
fi


bgzip "TMP/${params.prefix?:""}${genomeId}.spliceAI.scores.tsv"
tabix --force -s 1 -b 2 -e 2 "TMP/${params.prefix?:""}${genomeId}.spliceAI.scores.tsv.gz"

bgzip "TMP/${params.prefix?:""}${genomeId}.spliceAI.scores.merged.tsv"
tabix --force -s 1 -b 2 -e 2 "TMP/${params.prefix?:""}${genomeId}.spliceAI.scores.merged.tsv.gz"

mv TMP/*.gz ./
mv TMP/*.tbi ./
"""
}
