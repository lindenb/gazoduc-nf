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
include {moduleLoad} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'

String normChr(c) {
	if(c.startsWith("chr")) return c.substring(3);
	return c;
	}


workflow BEAGLE01 {
take:
	meta
	genomeId
	rows
main:
	version_ch = Channel.empty()

	jar_ch = DOWNLOAD_BEAGLE_JAR([:])
	version_ch = version_ch.mix(jar_ch.version)

	gm_ch = DOWNLOAD_GENETIC_MAPS([:], genomeId)
	version_ch = version_ch.mix(gm_ch.version)


	rows2_ch = rows.combine( gm_ch.output.splitCsv(sep:'\t',header:false) ).
		filter{T->T[0].contig.equals(T[1])}.
		map{T->T[0].plus(gmap:T[2])}
	

	if(params.genomes[genomeId].ucsc_name.equals("hg19")) {
		contigs1_ch =  Channel.from("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X" /* no chrY */) 
		contigs2_ch = rows.map{T->T.contig}.unique()

		contigs3_ch = contigs1_ch.combine(contigs2_ch).
			filter{T->normChr(T[0]).equals(normChr(T[1]))}.
			map{T->T[0]}.
			unique()

		ch2_ch = DOWNLOAD_REF_HG19_VCF([:], genomeId, contigs3_ch )
		version_ch = version_ch.mix(ch2_ch.version)

		rows3_ch = rows2_ch.combine( ch2_ch.output.splitCsv(sep:'\t',header:false) ).
			filter{T->T[0].contig.equals(T[1])}.
			map{T->T[0].plus(ref:T[2])}

		}
	else if(params.genomes[genomeId].ucsc_name.equals("hg38")) {
		ch1_ch = DOWNLOAD_REF_HG38_MANIFEST([:])
		version_ch = version_ch.mix(ch1_ch.version)

		ch2_ch = DOWNLOAD_REF_HG38_VCF([:], genomeId, ch1_ch.output.splitCsv(sep:'\t',header:false) )
		version_ch = version_ch.mix(ch2_ch.version)

		rows3_ch = rows2_ch.combine( ch2_ch.output.splitCsv(sep:'\t',header:false) ).
			filter{T->T[0].contig.equals(T[1][0])}.
			map{T->T[0].plus(ref:T[1][1])}
		}
	else
		{
		throw new IllegalArgumentException("undefined build for beagle")
		}

	bgl_ch = APPLY_BEAGLE([:], jar_ch.jar, genomeId, rows3_ch)	
	version_ch = version_ch.mix(bgl_ch.version)
	
	out_ch = bgl_ch.output.map{T->T[0].plus(phased_vcf:T[1])}

	version_ch = MERGE_VERSION("beagle",version_ch.collect())

emit:
	version = version_ch
	output = out_ch
}




process DOWNLOAD_BEAGLE_JAR {
tag "${params.beagle_jar_url}"
input:
	val(meta)
output:
	path("beagle.jar"),emit:jar
	path("version.xml"),emit:version
script:
def url = params.beagle_jar_url
"""
hostname 1>&2

wget -O beagle.jar "${url}"

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description"></entry>
</properties>
EOF
"""
}


process DOWNLOAD_REF_HG38_MANIFEST {
input:
	val(meta)
output:
	path("urls.tsv"),emit:output
	path("version.xml"),emit:version
script:
def base="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"
"""
wget -O - "${base}/phased-manifest_July2021.tsv" |\
	cut -f 1 |\
	grep '.vcf.gz\$' |\
	awk '{printf("${base}\t%s\\n",\$0);}' > urls.tsv
"""
}

process DOWNLOAD_REF_HG19_VCF {
tag "chr${contig}"
afterScript "rm -rf TMP"
maxForks 5
cpus 5
memory "2g"
input:
	val(meta)
	val(genomeId)
	val(contig)
output:
	path("chrom2vcf.txt"),emit:output
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
	def base= "https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf"
	def vcf = "chr${contig}.1kg.phase3.v5a.vcf.gz"
"""
hostname 1>&2
${moduleLoad("bcftools jvarkit")}
mkdir -p TMP
set -x

wget -O TMP/jeter.vcf.gz "${base}/${vcf}"

bcftools query  -f '%CHROM\\n'  TMP/jeter.vcf.gz | head -n 1 > TMP/jeter.a
cat TMP/jeter.a | java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr  --convert  RETURN_ORIGINAL -R "${reference}" > TMP/jeter.b

if cmp --silent TMP/jeter.a TMP/jeter.b ; then

	# no need to convert the chromosomes
	mv -v  TMP/jeter.vcf.gz "${vcf}"

else


	paste TMP/jeter.a TMP/jeter.b > TMP/jeter.c
	# rename chr wants index
	bcftools index --threads ${task.cpus}  TMP/jeter.vcf.gz
	bcftools annotate --threads ${task.cpus} -O z -o "${vcf}" --rename-chrs TMP/jeter.c  TMP/jeter.vcf.gz

fi

bcftools index --threads ${task.cpus} -t "${vcf}"
bcftools index -s "${vcf}" | cut -f 1 | tr "\\n" "\\t" > chrom2vcf.txt
echo "\${PWD}/${vcf}" >> chrom2vcf.txt

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description"></entry>
</properties>
EOF
"""
}


process DOWNLOAD_REF_HG38_VCF {
tag "${vcf}"
maxForks 5
input:
	val(meta)
	val(genomeId)
	tuple val(base),val(vcf)
output:
	path("chrom2vcf.txt"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}
wget -O "${vcf}" "${base}/${vcf}"
wget -O "${vcf}.tbi" "${base}/${vcf}.tbi"

bcftools index -s "${vcf}" | tr "\\n" "\\t" > chrom2vcf.txt
echo "\${PWD}/${vcf}" >> chrom2vcf.txt

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description"></entry>
</properties>
EOF
"""
}

workflow DOWNLOAD_GENETIC_MAPS {
	take:
		meta
		genomeId
	main:
		version_ch = Channel.empty()
		ch1 = GENETIC_MAPS_DOWNLOAD_ZIP([:],genomeId)
		version_ch = version_ch.mix(ch1.version)

		ch2 = GENETIC_MAPS_RENAME([:], genomeId , ch1.output.splitText().map{it.trim()} ) 
		version_ch = version_ch.mix(ch2.version)


		ch3 = GENETIC_MAPS_MERGE([:], ch2.output.collect())
		version_ch = version_ch.mix(ch3.version)
	emit:
		version = version_ch
		output = ch3.output
}



process GENETIC_MAPS_DOWNLOAD_ZIP {
tag "${genomeId}"
afterScript "rm -rf TMP OUT/plink.README.txt OUT/README.txt"
memory "2g"
input:
	val(meta)
	val(genomeId)
output:
	path("maps.list"),emit:output
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def url = genome.beagle_genetic_map
"""
hostname 1>&2
mkdir -p TMP

wget -O TMP/jeter.zip "${url}"

(cd TMP && unzip jeter.zip)

mv -v TMP OUT

find \${PWD}/OUT -type f -name "plink*.map" > maps.list

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description"></entry>
</properties>
EOF
"""
}

process GENETIC_MAPS_RENAME {
tag "${mapFile}"
afterScript "rm -rf TMP"
memory "2g"
input:
        val(meta)
        val(genomeId)
	val(mapFile)
output:
       	path("chrom2map.tsv"),emit:output
        path("version.xml"),emit:version
script:
        def genome = params.genomes[genomeId]
        def reference = genome.fasta
	def fname = file(mapFile).name
"""
hostname 1>&2
${moduleLoad("jvarkit")}
mkdir -p TMP

java  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${reference}" --column 1 --convert SKIP '${mapFile}' > TMP/${fname}
mv  TMP/${fname} ./

echo "\${PWD}/${fname}"	| awk -F '.' '{printf("%s\t%s\\n",\$(NF-2),\$0);}' |\
	java  -Xmx${task.memory.giga}g -Djava.io.tmpdir=.  -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${reference}" --column 1 --convert SKIP > chrom2map.tsv


##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description"></entry>
</properties>
EOF
"""
}

process GENETIC_MAPS_MERGE {
executor "local"
tag "N=${L.size()}"
input:
	val(meta)
	val(L)
output:
       	path("chrom2map.tsv"),emit:output
        path("version.xml"),emit:version
"""
cat ${L.join(" ")} | sort -t '\t' -k1,1V  > chrom2map.tsv
test -s chrom2map.tsv
##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
       	<entry key="description"></entry>
</properties>
EOF
"""
}


process APPLY_BEAGLE {
tag "${row.interval}"
memory "5g"
cpus 5
afterScript "rm -rf TMP"
input:
	val(meta)
	path(beagle)
	val(genomeId)
	val(row)
output:
	tuple val(row),path("phased.bcf"),emit:output
	path("version.xml"),emit:version
script:
	if(!row.containsKey("vcf")) throw new IllegalArgumentException("row.vcf is missing")
	if(!row.containsKey("interval")) throw new IllegalArgumentException("row.interval is missing")
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP

if ${row.vcf.endsWith(".bcf")} ; then
	bcftools view  --threads ${task.cpus}  --regions "${row.interval}" "${row.vcf}" -O z -o TMP/jeter.vcf.gz
	bcftools index  --threads ${task.cpus}  -f -t TMP/jeter.vcf.gz
fi


java  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar ${beagle} \
	ref=${row.ref} \
	map=${row.gmap} \
	chrom=${row.interval} \
	gt=${row.vcf.endsWith(".bcf")?"TMP/jeter.vcf.gz":"${row.vcf}"} \
	out=TMP/phased \
	nthreads=${task.cpus} 1>&2

# add dictionary
bcftools reheader --fai "${reference}.fai" -T TMP/tmp --threads ${task.cpus} -o TMP/phased2.vcf.gz  TMP/phased.vcf.gz 
mv -v TMP/phased2.vcf.gz TMP/phased.vcf.gz

bcftools view -O b -o phased.bcf TMP/phased.vcf.gz
bcftools index phased.bcf

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description"></entry>
</properties>
EOF
"""
}
