
process CARDIOPANEL_DOWNLOAD {
tag "${meta1?:fasta.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta4),path(gtf),path(gtf_tbi)
output:
    tuple val(meta1),path("*IN.bed.gz"), path("*IN.bed.gz.tbi"), path("*IN.header"), emit:bed
    tuple val(meta1),path("*OUT.bed.gz"),path("*OUT.bed.gz.tbi"),path("*OUT.header"),emit:bed_extended
    path("versions.yml"),emit:versions
script:
    def TAG = task.ext.tag?:"CARDIOPANEL"
    def FILE = task.ext.file?:"\${HOME}/notebook/data/genes/20220421.genes.cardioPanel_CECAS_CHU_Nantes.txt"
    def WHATIZ = "Cardiopanel ${FILE}"
	def extend = task.ext.extend?:1000
"""
hostname 1>&2
mkdir -p TMP
set -o pipefail
mkdir -p TMP

export LC_ALL=C

if test -f "${FILE}"
then
	grep -v "#" "${FILE}" |\\
		grep -v '^\$' |\\
		sort -T TMP |  uniq > TMP/genes.txt


	echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="${WHATIZ}">' > ${TAG}_IN.header
	echo '##INFO=<ID=${TAG}_NEAR,Number=.,Type=String,Description="Near gene distance=${extend}. ${WHATIZ}">' > ${TAG}_OUT.header

else

	touch TMP/genes.txt
	
	echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="${TAG} NOT AVAILABLE">' > ${TAG}_IN.header
	echo '##INFO=<ID=${TAG}_NEAR,Number=.,Type=String,Description="NOT AVAILABLE">' > ${TAG}_OUT.header

fi

jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP gtf2bed  \\
		--columns "gtf.feature,gene_name" -R "${fasta}"  "${gtf}" |\\
	awk -F '\t' '\$4=="gene" && \$5!="." && \$5!=""' |\\
	grep -f  TMP/genes.txt -F |\\
	cut -f1,2,3,5 |\\
    LC_ALL=C sort --buffer-size=${task.memory.mega}M -t '\t' -k4,4 -T TMP  |\\
	uniq > TMP/genes.bed


join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,1.4' TMP/genes.bed TMP/genes.txt |\\
	sort -T TMP -t '\t' -k1,1 -k2,2n |\\
	uniq |\\
	bgzip > TMP/${TAG}_IN.bed.gz

tabix --force -p bed TMP/${TAG}_IN.bed.gz


gunzip -c TMP/${TAG}_IN.bed.gz |\\
	bedtools slop -b ${extend} -g ${fai} |\\
	LC_ALL=C sort -t '\t' -k1,1 -k2,2n -T TMP > TMP/extended.bed


bedtools subtract -a TMP/extended.bed -b  TMP/${TAG}_IN.bed.gz |\\
	awk -F '\t' 'int(\$2) < int(\$3)' |\\
	LC_ALL=C sort -t '\t' -k1,1 -k2,2n -T TMP |\\
	uniq |\\
	bgzip  > TMP/${TAG}_OUT.bed.gz

tabix --force -p bed TMP/${TAG}_OUT.bed.gz

mv TMP/*.bed.gz ./
mv TMP/*.bed.gz.tbi ./




cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${FILE}"
END_VERSIONS
"""
}
