
process GREENDB_DOWNLOAD {
tag "${meta1.id?:fasta.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
output:
	tuple val(meta1),path("*.bed.gz"),path("*.bed.gz.tbi"),path("*.header"),emit:bed
    path("versions.yml"),emit:versions
script:
    def TAG = "GREENDB"
    if(!meta1.ucsc_name) throw new IllegalArgumentException("${task.process} undefined ucsc_name");
    def url="";
    if(meta1.ucsc_name.equals("hg38")) {
        url = "https://zenodo.org/record/5636209/files/GRCh38_GREEN-DB.bed.gz?download=1";
    } else if(meta1.ucsc_name.equals("hg19")) {
        url = "https://zenodo.org/record/5636209/files/GRCh37_GREEN-DB.bed.gz?download=1";
    } else {
        throw new IllegalArgumentException("${task.process} unknown ucsc_name");
        }

    def whatis = "GREEN-DB is a comprehensive collection of 2.4 million regulatory elements in the human genome collected from previously published databases, high-throughput screenings and functional studies. ${url}"

"""

hostname 1>&2
mkdir -p TMP
env | grep -i proxy 1>&2


curl -o TMP/jeter.bed.gz -L "${url}"

gunzip -c TMP/jeter.bed.gz |\\
	grep -v '^chromosome' |\\
	cut -f1-3,5 |\\
	jvarkit bedrenamechr -f "${fasta}" --column 1 --convert SKIP  |\\
	LC_ALL=C sort  --buffer-size=${task.memory.mega}M  -T . -t '\t' -k1,1 -k2,2n |\\
	uniq |\\
	bgzip > TMP/${TAG}.bed.gz


tabix --force -p bed TMP/${TAG}.bed.gz

mv TMP/${TAG}.bed.gz ./
mv TMP/${TAG}.bed.gz.tbi ./

echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="${whatis}. ${url}.">' > ${TAG}.header


cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${url}"
END_VERSIONS
"""
}
