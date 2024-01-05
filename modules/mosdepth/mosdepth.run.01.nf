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
include {parseBoolean;isBlank} from '../utils/functions.nf'

process MOSDEPTH_RUN_01 {
	tag "${row.sample} ${file(row.bam).name}"
	afterScript "rm -rf TMP"
	/** cpus 1 must be specified in config */
	input:
		val(meta)
		val(executable)
		val(row)
	output:
		tuple val(row),path("${row.sample}${row.suffix?:""}.mosdepth.global.dist.txt"),emit:globaldist
		tuple val(row),path("${row.sample}${row.suffix?:""}.mosdepth.region.dist.txt"),optional:true,emit:regiondist
		tuple val(row),path("${row.sample}${row.suffix?:""}.mosdepth.summary.txt"),emit:summary
		tuple val(row),path("${row.sample}${row.suffix?:""}.per-base.bed.gz"),optional:true,emit:perbase
		tuple val(row),path("${row.sample}${row.suffix?:""}.regions.bed.gz"),optional:true,emit:regions
		tuple val(row),path("${row.sample}${row.suffix?:""}.output.tsv"),emit:output
		path("version.xml"),emit:version
	script:
		if(!row.containsKey("sample")) throw new IllegalArgumentException("[mosdepth]row.sample is missing");
		if(!row.containsKey("reference") && !row.containsKey("genomeId")) throw new IllegalArgumentException("[mosdepth]row.reference/genomeId is missing");

		def reference

		if( row.containsKey("reference") ) {
			reference = row.reference 
			}
		else 	
			{
			reference = params.genomes[row.genomeId].fasta
			}
		def bed = row.bed?row.bed:file("NO_FILE")
		def mapq = row.mapq?:params.mosdepth.mapq
		def suffix = (row.sample?:"")+(row.suffix?:"")
		def extra = params.mosdepth.args
	"""
	hostname 1>&2
	mkdir -p TMP

        ${executable} ${bed.name.equals("NO_FILE")?"":"--by \"${bed}\""} \
 		-t ${task.cpus} --fasta "${reference}" --mapq ${mapq} ${extra} \
		'TMP/${suffix}' "${row.bam}"

	mv -v TMP/${suffix}* ./

cat << EOF | paste -s > output.tsv
sample
reference
bam
globaldist
regiondist
summary
perbase
regions
EOF


echo -ne "${row.sample}\t${reference}\t${row.bam}\t\${PWD}/${suffix}.mosdepth.global.dist.txt" >> output.tsv

if  test -f "${suffix}.mosdepth.region.dist.txt" ; then
	echo -ne "\t\${PWD}/${suffix}.mosdepth.region.dist.txt" >> output.tsv
else
	echo -ne "\tNO_FILE" >> output.tsv
fi

echo -en "\t\${PWD}/${suffix}.mosdepth.summary.txt" >> output.tsv


if  test -f "${suffix}.per-base.bed.gz" ; then
	echo -ne "\t\${PWD}/${suffix}.per-base.bed.gz" >> output.tsv
else
	echo -ne "\tNO_FILE" >> output.tsv
fi

if  test -f "${suffix}.regions.bed.gz" ; then
	echo -ne "\t\${PWD}/${suffix}.regions.bed.gz" >> output.tsv
else
	echo -ne "\tNO_FILE" >> output.tsv
fi

echo >> output.tsv

mv output.tsv "${suffix}.output.tsv"

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">apply mosdepth to bam</entry>
	<entry key="bam">${row.bam}</entry>
        <entry key="mosdepth.version">\$(${executable} --version)</entry>
        <entry key="mapq">${mapq}</entry>
</properties>
EOF
	"""
	}
