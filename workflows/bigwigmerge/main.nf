
workflow {
	MERGE_BIGWIG(params.fasta, file(params.samplesheet))
	}


workflow MERGE_BIGWIG {
	take:
		fasta
		samplesheet
	main:
		fai = file(fasta+".fai")
		dict = file("${file(fasta).getParent()}/${file(fasta).getSimpleName()}.dict")



		input_ch = Channel.fromPath(samplesheet).
			splitCsv(header:true,sep:',').
			map{T->{
			if(T.containsKey("track")) return T;
			return T.plus("track":params.default_track_name);
			}}.
			map{T->[T.track,T.bigWig]}.
			groupTuple()

		contigs_ch = Channel.fromPath(fai).
			splitCsv(header:false,sep:'\t').
			map{[it[0],it[1]]}

		exec_ch = GET_BEDGRAPH_TO_BIGWIG()
		
		ch1_ch = contigs_ch.combine(input_ch)
		merge1_ch = MERGE1(fasta,fai,dict,ch1_ch)
		merge1_ch = MERGE2(fasta,fai,dict,exec_ch.output,merge1_ch.output.
			map{[it[0],it[1]+"\t"+it[2]]}.
			groupTuple())

	}

process MERGE1 {
tag "${contig} ${track} N=${L.size()}"
memory "3g"
afterScript "rm -rf TMP"
input:
	path(fasta)
	path(fai)
	path(dict)
	tuple val(contig),val(contigLen),val(track),val(L)
output:
	tuple val(track),val(contig),path("${track}.${contig}.bed.gz"),emit:output
"""
hostname 1>&2
mkdir -p TMP
set -o pipefail
module load jvarkit

cat << EOF > TMP/jeter.list
${L.join("\n")}
EOF

java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${HOME}/jeter8.jar -R '${fasta}' --interval "${contig}:1-${contigLen}" TMP/jeter.list -m ${params.method} | gzip --best > TMP/jeter.bed.gz
mv -v TMP/jeter.bed.gz ${track}.${contig}.bed.gz
"""
}


process GET_BEDGRAPH_TO_BIGWIG {
output:
	
	path("bedGraphToBigWig"),emit:output
"""
wget -O bedGraphToBigWig "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bedGraphToBigWig"
chmod +x bedGraphToBigWig
"""
}

process MERGE2 {
tag "${track} N=${L.size()}"
memory "3g"
afterScript "rm -rf TMP"
input:
	path(fasta)
	path(fai)
	path(dict)
	path(bedGraphToBigWig)
	tuple val(track),val(L)
output:
	path("${params.prefix}${track}.bigWig"),emit:output
"""
hostname 1>&2
mkdir -p TMP
set -o pipefail

cat << EOF | sort -T TMP -k1,1 -t '\t' | cut -f 2 >  TMP/jeter.list
${L.join("\n")}
EOF

xargs -a TMP/jeter.list gunzip -c > TMP/jeter.bed
./${bedGraphToBigWig} TMP/jeter.bed ${fai}  TMP/jeter.bw
mv -v TMP/jeter.bw "${params.prefix}${track}.bigWig"
"""
}

