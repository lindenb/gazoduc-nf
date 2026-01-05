/*

Copyright (c) 2026 Pierre Lindenbaum

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

include {moduleLoad;getVersionCmd} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'


workflow DOWNLOAD_GNOMAD_SV_01 {
     take:
        meta /* meta */
        genomeId /* indexed fasta reference either hg38 of hg39*/
     main:
        version_ch = Channel.empty()


	if(params.genomes[genomeId].ucsc_name.equals("hg19")) {
		ch0 = DOWNLOAD_BED_HG19(meta)
		version_ch = version_ch.mix(ch0.version)

		ch1 = RENAME_CONTIGS_HG19([:], ch0.output, genomeId)
		version_ch = version_ch.mix(ch1.version)

		dest_bed = ch1.bed
		dest_idx = ch1.index
		}
	else if(params.genomes[genomeId].ucsc_name.equals("hg38")) {

		ch0 = DOWNLOAD_BED_HG38([:],genomeId)
		version_ch = version_ch.mix(ch0.version)


		dest_bed = ch0.bed
		dest_idx = ch0.index
		}
	else
		{
		ch3 = OTHER_REF([:], ch0.output)

		version_ch = version_ch.mix(ch3.version)

		dest_bed = ch3.bed
		dest_idx = ch3.index
		}
		
	version_ch = MERGE_VERSION("GnomadSV",version_ch.collect())

	emit:
		bed = dest_bed
		index = dest_idx
		version = version_ch
	}

process DOWNLOAD_BED_HG19 {
input:
	val(meta)
output:
	path("gnomad.sv.37.bed.gz"),emit:output
	path("version.xml"),emit:version
script:
	
	def genome = params.genomes["hs37d5"]
	def url = genome.gnomad_sv_url
"""
wget -O "gnomad.sv.37.bed.gz" "${url}" 

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">Download gnomad sv</entry>
        <entry key="gnomad.sv.url"><a>${url}</a></entry>
        <entry key="versions">${getVersionCmd("wget")}</entry>
</properties>
EOF
"""
}


process DOWNLOAD_BED_HG38 {
tag "${genomeId}"
afterScript "rm -rf TMP"
memory "5g"
input:
	val(meta)
	val(genomeId)
output:
	path("gnomad.38.sv.bed.gz"),emit:bed
	path("gnomad.38.sv.bed.gz.tbi"),emit:index
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
"""
hostname 1>&2
${moduleLoad("bcftools htslib jvarkit")}
set -x


mkdir -p TMP
for C in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
	wget -q -O "TMP/chr\${C}.vcf.gz"  "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/genome_sv/gnomad.v4.0.sv.chr\${C}.vcf.gz"
	bcftools index -f "TMP/chr\${C}.vcf.gz"
 
	echo "TMP/chr\${C}.vcf.gz" >> TMP/jeter.list
done

test -s TMP/jeter.list



bcftools concat --file-list TMP/jeter.list -O u | bcftools view --header-only -O b -o TMP/header.bcf

echo -n "%CHROM\t%POS0\t%END\t%REF\t%ALT\t%FILTER" > TMP/jeter.query

bcftools view --header-only TMP/header.bcf | grep "^##INFO=" | tr "<>" "\\n" | grep "^ID=" | sed 's/,.*//' | cut -d '=' -f 2- | awk '{printf("\t%INFO/%s",\$1);} END {printf("\\n");}' >> TMP/jeter.query

sed 's%INFO/%%g' < TMP/jeter.query | tr "a-z" "A-Z" | sed 's/CHROM\tPOS0\tEND/chrom\tstart\tend/'  | tr -d '%' | awk '{printf("#%s\\n",\$0);}' > TMP/gnomad.38.sv.bed

QUERY=`cat TMP/jeter.query`

echo "\${QUERY}" 1>&2

cat TMP/jeter.list | while read F
do
	bcftools query -f "\${QUERY}"  "\${F}" |\\
		java -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${reference}" --column 1 --convert SKIP >> TMP/gnomad.38.sv.bed
done

test -s TMP/gnomad.38.sv.bed

bgzip -f TMP/gnomad.38.sv.bed
tabix --comment '#'  -p bed -f TMP/gnomad.38.sv.bed.gz

mv TMP/gnomad.38.sv.bed.gz ./
mv TMP/gnomad.38.sv.bed.gz.tbi ./

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">download for hg38</entry>
</properties>
EOF
"""
}

process RENAME_CONTIGS_HG19 {
tag "${genomeId}"
afterScript "rm -rf TMP"
input:
	val(meta)
	path(bed)
	val(genomeId)
output:
	path("gnomad.37.sv.bed.gz"),emit:bed
	path("gnomad.37.sv.bed.gz.tbi"),emit:index
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
"""
${moduleLoad("jvarkit htslib")}
## set -o pipefail non gunzip...
mkdir -p TMP

gunzip -c "${bed}" | head -n1  > TMP/gnomad.37.sv.bed


gunzip -c "${bed}" |\
	tail -n+2 |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n >> TMP/gnomad.37.sv.bed

bgzip -f TMP/gnomad.37.sv.bed
tabix --comment '#'  -p bed -f TMP/gnomad.37.sv.bed.gz

mv TMP/gnomad.37.sv.bed.gz ./
mv TMP/gnomad.37.sv.bed.gz.tbi ./

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">fix chromosomes for hg19</entry>
</properties>
EOF
"""
}



process OTHER_REF {
input:
	val(meta)
	path(bed)
output:
	path("gnomad.x.sv.bed.gz"),emit:bed
	path("gnomad.x.sv.bed.gz.tbi"),emit:index
	path("version.xml"),emit:version
script:
	def url2 = meta.gnomad_sv_grch37_bed_url
"""
${moduleLoad("htslib")}
gunzip -c "${bed}" | head -n 1 > gnomad.x.sv.bed
bgzip -f gnomad.x.sv.bed
tabix --comment '#'  -p bed -f  gnomad.x.sv.bed.gz


#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">Download gnomad sv but undefined ref</entry>
</properties>
EOF
"""
}
