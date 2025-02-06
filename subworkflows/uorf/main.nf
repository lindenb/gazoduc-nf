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
include {getKeyValue;getModules} from '../../modules/utils/functions.nf'
include {k1_signature} from '../../modules/utils/k1.nf'
//include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {VCF_TO_BED} from '../vcf2bed'


def k1 = k1_signature()

workflow UORF {
	take:
		reference
		vcf
		userbed
		samples
	main:
		version_ch  = Channel.empty()



		gencode_ch = DOWNLOAD_GENCODE(reference)
		kg2bed_ch = KG2BED(gencode_ch.output)
		kg2gff_ch = KG2GFF3(gencode_ch.output)
		uorfgff_ch = UORFGFF3(reference, gencode_ch.output)

		vcf2bed_ch = VCF_TO_BED(vcf)
		vcf2bed_ch.bed.splitCsv(header:false, sep:'\t').
			map{[it[0]+":"+((it[1] as int)+1)+"-"+it[2],file(it[3]),file(it[4])]}.
			set{rgn_vcf}

		annot_ch = ANNOTATE(reference, userbed, kg2bed_ch.output, kg2gff_ch.output, uorfgff_ch.output, samples, rgn_vcf)

		merge_ch = MERGE(annot_ch.output.collect())
	emit:
		output = merge_ch.output
	}


process DOWNLOAD_GENCODE {
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
        path(genome)
output:
        path("wgEncodeGencode.txt.gz"),emit:output
script:
	def fai = genome.find{it.name.endsWith(".fai")}
	def dict = genome.find{it.name.endsWith(".dict")}
"""
hostname 1>&2
mkdir -p TMP
set -o pipefail

cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg38}\thttps://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeBasicV47.txt.gz
1:${k1.hg19}\thttps://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeBasicV19.txt.gz
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sed 's/^chr//' | sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv

join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv | sort | uniq > TMP/jeter.url

echo -n "URL is " 1>&2 && cat TMP/jeter.url 1>&2
test -s TMP/jeter.url

wget -O - `cat TMP/jeter.url` |\\
	gunzip -c |\\
	awk -F '\t' '(\$7!=\$8)' |\\
	jvarkit  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP bedrenamechr -R "${dict}" --column 3 --convert SKIP  |\\
	gzip > TMP/wgEncodeGencode.txt.gz
mv TMP/wgEncodeGencode.txt.gz ./
"""
}

process KG2BED {
tag "${kg.name}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	path(kg)
output:
	path("${kg.baseName}.bed"),emit:output
script:
"""
mkdir -p TMP
set -o pipefail

gunzip -c "${kg}" |\\
	jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP kg2bed |\\
	grep -F -w UTR5 |\\
	cut -f1,2,3 |\\
	sort -T TMP -k1,1 -k2,2n -t '\t' |\\
	bedtools merge > TMP/jeter.bed

mv TMP/jeter.bed "${kg.baseName}.bed"
"""
}


process KG2GFF3 {
tag "${kg.name}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	path(kg)
output:
	path("${kg.baseName}.gff3.gz"),emit:output
script:
"""
mkdir -p TMP
set -o pipefail

gunzip -c "${kg}" |\\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${HOME}/jvarkit.jar kg2gff |\\
	LC_ALL=C sort -T TMP -k1,1 -k4,4n -t '\t' |\\
	gzip > TMP/jeter.gff3.gz

mv TMP/jeter.gff3.gz ./${kg.baseName}.gff3.gz
"""
}


process UORFGFF3 {
tag "${kg.name}"
label "process_medium"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	path(genome)
	path(kg)
output:
	path("${kg.baseName}.uorf.gff3.gz"),emit:output
script:
	def fasta = genome.find{it.name.endsWith("a")}
"""
mkdir -p TMP
set -o pipefail

gunzip -c "${kg}" |\\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${HOME}/jvarkit.jar gff3upstreamorf -R '${fasta}' |\\
	LC_ALL=C sort -T TMP -k1,1 -k4,4n -t '\t' |\\
	gzip > TMP/jeter.gff3.gz

mv TMP/jeter.gff3.gz ./${kg.baseName}.uorf.gff3.gz
"""
}

process ANNOTATE {
tag "${rgn} ${vcf.name}"
label "process_medium"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	path(reference)
	path(userbed)
	path(bed)
	path(gff3)
	path(uorfgff3)
	path(samples)
	tuple val(rgn), path(vcf),path(vcfix)
output:
	path("uorf.*"),emit:output
script:
	def acn = "SO:0001629,SO:0001818"
	def fasta = reference.find{it.name.endsWith("a")}
	def fai = reference.find{it.name.endsWith("fai")}
"""
mkdir -p TMP
set -o pipefail
set -x

echo '${rgn}' | awk -F '[:-]' '{printf("%s\t%d\t%s\\n",\$1,int(\$2)-1,\$3);}' > TMP/jeter.a

bedtools intersect -a TMP/jeter.a -b '${bed}' |\\
	sort -T TMP -k1,1 -k2,2n -t '\t' > TMP/jeter.bed

echo -e 'INFO/BCSQ\tSTD_BCSQ' > TMP/rename_annot.txt

if ${userbed.name.contains(".")}
then
	bedtools intersect -a '${userbed}' -b TMP/jeter.bed |\\
		sort -T TMP -k1,1 -k2,2n -t '\t' > TMP/jeter2.bed
	mv TMP/jeter2.bed TMP/jeter.bed
fi


if test ! -s TMP/jeter.bed
then
	tail -1 "${fai}" | awk -F '\t' '{printf("%s\t0\t1\\n",\$1);}' > TMP/jeter.bed
fi

if ${samples.name.contains(".")}
then
	bcftools query -l '${vcf}' | sort | uniq > TMP/sn.a
	sort -T TMP -t '\t' -k1,1 '${samples}' | uniq > TMP/sn.b
	join -t '\t' -1 1 -2 1 -o '2.1,2.2' TMP/sn.a TMP/sn.b > TMP/samples_pheno.tsv
	cut -f 1 TMP/samples_pheno.tsv > TMP/samples.txt
	awk -F '\t' '(\$2=="controls")' TMP/samples_pheno.tsv | cut -f1 | sort | uniq > TMP/controls.txt
	awk -F '\t' '(\$2=="case")' TMP/samples_pheno.tsv | cut -f1 | sort | uniq > TMP/cases.txt
fi

bcftools csq \\
	${samples.name.contains(".")?"--samples-file TMP/samples.txt":""} \\
	--regions-file TMP/jeter.bed \\
	--local-csq \\
	--ncsq 1000 \\
	--gff-annot "${gff3}" \\
	--fasta-ref "${fasta}" \\
	-O v \\
	"${vcf}" |\\
	jvarkit -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP vcffilterso --acn "${acn}" --filterin "BASIC_SO" |\\
	bcftools view -O u -e 'FILTER ~ "BASIC_SO"' |\\
	bcftools annotate -x 'FILTER/BASIC_SO' -O b -o TMP/jeter.bcf

# remove previous annot rename-annot is buggy with current bcftools ?
# bcftools annotate --rename-annots TMP/rename_annot.txt -O b -o TMP/jeter2.bcf TMP/jeter.bcf
bcftools annotate -x 'INFO/BCSQ' -O b -o TMP/jeter2.bcf TMP/jeter.bcf
mv  TMP/jeter2.bcf TMP/jeter.bcf

bcftools index -f TMP/jeter.bcf

bcftools csq \\
	--regions-file TMP/jeter.bed \\
	--local-csq \\
	--ncsq 1000 \\
	--gff-annot "${uorfgff3}" \\
	--fasta-ref "${fasta}" \\
	TMP/jeter.bcf |\
	jvarkit -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP vcffilterso --acn "${acn}" |\\
	bcftools view -O b -o TMP/jeter2.bcf

if test -s TMP/cases.txt && test -s TMP/controls.txt
then
	bcftools +contrast -0 TMP/controls.txt -1 TMP/cases.txt -a PASSOC,FASSOC,NASSOC,NOVELAL,NOVELGT -O b -o TMP/jeter3.bcf TMP/jeter2.bcf
	mv TMP/jeter3.bcf TMP/jeter2.bcf
fi

bcftools index -f TMP/jeter2.bcf

MD5=`md5sum TMP/jeter2.bcf | awk '{print \$1}'`

mv TMP/jeter2.bcf "uorf.\${MD5}.uorf.bcf"
mv TMP/jeter2.bcf.csi "uorf.\${MD5}.uorf.bcf.csi"

"""
}



process MERGE {
label "process_medium"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	path("VCF/*")
output:
	tuple path("concat.bcf"),path("concat.bcf.csi"),emit:output
script:
"""
mkdir -p TMP
find VCF -type l -name "*.bcf" > TMP/jeter.list

bcftools concat --threads "${task.cpus}" -a --file-list TMP/jeter.list -O b -o TMP/concat.bcf
bcftools index --threads "${task.cpus}" -f TMP/concat.bcf
mv TMP/concat.bcf ./
mv TMP/concat.bcf.csi ./
"""
}

