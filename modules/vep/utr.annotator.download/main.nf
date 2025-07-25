process DOWNLOAD_UTR_ANNOTATOR {
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
output:
    tuple val(meta1),path("*.txt"),emit:output
    path("versions.yml"),emit:versions
script:
    def f1=""
    if((!meta1.ucsc_name)) {
        throw new IllegalArgumentException("${task.process} missing ucsc_name");
    } else if(meta1.ucsc_name.equals("hg38")) {
        f1="uORF_5UTR_GRCh38_PUBLIC.txt"
    } else if(meta1.ucsc_name.equals("hg19")) {
        f1="uORF_5UTR_GRCh37_PUBLIC.txt"
    } else {
        throw new IllegalArgumentException("${task.process} undefined ucsc_name");
    }
   def url = "https://github.com/ImperialCardioGenetics/UTRannotator/raw/refs/heads/master/${f1}"
"""
curl -L -o "${f1}" "${url}"

cat << END_VERSIONS > versions.yml
"${task.process}":
    url: "${url}"
END_VERSIONS
"""
}

