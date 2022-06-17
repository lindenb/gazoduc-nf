include { getKeyValue; getModules} from '../../modules/utils/functions.nf'

process SAMTOOLS_SAMPLES01 {
tag "${file(bams).name}"
afterScript "rm -f jeter.txt jeter.tsv"
input:
	val(meta)
	val(reference)
	val(bams)
output:
	path("sample2bam.tsv"),emit:out
	path("version.xml"),emit:version
script:
"""
hostname 2>&1
module load ${getModules("samtools")}
set -o pipefail


samtools samples -f "${reference}" < ${bams} | sort | uniq > jeter.tsv 

# no empty samples
awk -F '\t' '(\$1==".")' jeter.tsv > jeter.txt 
test ! -s jeter.txt

# no empty ref
awk -F '\t' '(\$3==".")' jeter.tsv > jeter.txt 
test ! -s jeter.txt

# no dup samples
cut -f 1 jeter.tsv | sort -T . | uniq -d > jeter.txt
test ! -s jeter.txt

# no dup bam
cut -f 2 jeter.tsv | sort -T . | uniq -d > jeter.txt
test ! -s jeter.txt


cut -f1,2 jeter.tsv | sort -T . -t '\t' -k1,1  > sample2bam.tsv

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">extract sample names from BAM metadata</entry>
	<entry key="input">${bams}</entry>
        <entry key="samtools">\$(samtools  --version | head -n 1| cut -d ' ' -f2)</entry>
</properties>
EOF
"""

}

