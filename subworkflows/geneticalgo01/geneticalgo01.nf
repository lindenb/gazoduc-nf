include {moduleLoad} from '../../modules/utils/functions.nf'

if(!params.containsKey("step_id")) throw new IllegalArgumentException("params.step_id missing");
if(!params.containsKey("step_name")) throw new IllegalArgumentException("params.step_name missing");

workflow GENETICALGO {
	take:
		genomeId
        	vcf
		cases
		controls
		exclude_bed_ch /** array  [selid, vcf] */
	main:
		intervals_ch = Channel.empty()		
		


		sel0a_ch = Channel.of(
			["exon10","Exons slop 10",10,"(\$3==\"exon\")"],
			["exon100","Exons slop 100",100,"(\$3==\"exon\")"],
			["gene10","Genes slop 10",10,"(\$3==\"gene\")"],
			["gene100","Genes slop 100",100,"(\$3==\"gene\")"],
			["gene500","Genes slop 500",500,"(\$3==\"gene\")"]
			).
			map{T->[ selid:T[0], selname:T[1], slop:T[2], selexpr:T[3] ]}

		sel0b_ch  = SELECT01(
			genomeId,
			sel0a_ch.map{T->[T.selid,T]}.
				join(exclude_bed_ch,remainder:true).
				filter{it[1]!=null}.
				map{T->T[1].plus("exclude_bed":T[2])}
			)

		intervals_ch = intervals_ch.mix(sel0b_ch.output)


		/** Whole genome */
		sel1a_ch = Channel.of([selid:"whole",selname:"All"])
	

		sel1b_ch  = WHOLE_GENOME(
			genomeId,
			sel1a_ch.map{T->[T.selid,T]}.
				join(exclude_bed_ch,remainder:true).
				filter{it[1]!=null}.
				map{T->T[1].plus("exclude_bed":T[2])}
			)
		
		intervals_ch = intervals_ch.mix(sel1b_ch.output)

		/** introns */
		sel2a_ch = Channel.of([selid:"introns",selname:"Introns"])
		sel2b_ch  = INTRONS(
			genomeId,
			sel2a_ch.map{T->[T.selid,T]}.
				join(exclude_bed_ch,remainder:true).
				filter{it[1]!=null}.
				map{T->T[1].plus("exclude_bed":T[2])}
			)

		intervals_ch = intervals_ch.mix(sel2b_ch.output)

		/** gene upstream/downstream */
		sel3a_ch = Channel.of(
			["upstream1kb","Upstream 1Kb","upstream",1_000],
			["upstream10kb","Upstream 10Kb","upstream",10_000],
			["upstream100kb","Upstream 100Kb","upstream",100_000],
                        ["downstream1kb","Downstream 1Kb","downstream",1_000],
                        ["downstream10kb","Downstream 10Kb","downstream",10_000],
                        ["downstream100kb","Downstream 100Kb","downstream",100_000],
			).
			map{T->[ selid:T[0], selname:T[1], side:T[2], extend:T[3]]}

		sel3b_ch = XXSTREAM (
			genomeId,
			sel3a_ch.map{T->[T.selid,T]}.
				join(exclude_bed_ch,remainder:true).
				filter{it[1]!=null}.
				map{T->T[1].plus("exclude_bed":T[2])}
			)
		intervals_ch = intervals_ch.mix(sel3b_ch.output)

		/** INTERGENIC **/

		sel4a_ch = Channel.of(0,10,100).
			map{T->[ selid:"intergenic"+T, selname:"Intergenic Gene slop "+T+"bp", slop:T]}

		sel4b_ch = INTERGENIC(
			genomeId,
			sel4a_ch.map{T->[T.selid,T]}.
				join(exclude_bed_ch,remainder:true).
				filter{it[1]!=null}.
				map{T->T[1].plus("exclude_bed":T[2])}
			)
		intervals_ch = intervals_ch.mix(sel4b_ch.output)

		/** regulation */
		reg_ch  = DOWNLOAD_GFF_REG(genomeId)
		sel5a_ch = reg_ch.types.splitText().
			map{it.trim()}.
			map{T->[ selid:T, selname:T]}

		sel5b_ch = REGULATORY(
			reg_ch.output,
			sel5a_ch.map{T->[T.selid,T]}.
				join(exclude_bed_ch,remainder:true).
				filter{it[1]!=null}.
				map{T->T[1].plus("exclude_bed":T[2])}
			)
		intervals_ch = intervals_ch.mix(sel5b_ch.output)


		/** first intron */
		sel6a_ch = Channel.of([selid:"firstintron",selname:"FirstIntron"])
		sel6b_ch  = FIRST_INTRON(
			genomeId,
			sel6a_ch.map{T->[T.selid,T]}.
				join(exclude_bed_ch,remainder:true).
				filter{it[1]!=null}.
				map{T->T[1].plus("exclude_bed":T[2])}
			)
		intervals_ch = intervals_ch.mix(sel6b_ch.output)

		apply_ch = APPLY(  genomeId,
			vcf,
			cases,
			controls,
			intervals_ch.map{T->T[0].plus(bed:T[1])}
			)

		to_mqc = apply_ch.mqc
		mqc_table_ch= MQC_TABLE(genomeId, apply_ch.tsv.map{T->T[0].selid+"\t"+T[1]}.collect())
		to_mqc = to_mqc.mix(mqc_table_ch.mqc)

	emit:
		mqc = to_mqc
		tsv = apply_ch.tsv.map{T->[T[0].plus(step_id:params.step_id,step_name:params.step_name),T[1]]}
		exclude_bed = apply_ch.exclude.map{T->[T[0].selid,T[1]]}
	}


process WHOLE_GENOME {
	tag "${row.selid}"
	afterScript "rm -rf TMP"
	input:
		val(genomeId)
		val(row)
	output:
		tuple val(row),path("select.bed"),emit:output
	script:
		def genome = params.genomes[genomeId]
		def fasta = genome.fasta
		def exclude_bed = row.exclude_bed==null?file("NO_FILE"):row.exclude_bed
	"""
	hostname 1>&2
	${moduleLoad("bedtools")}
	mkdir -p TMP
	awk -F '\t' '{printf("%s\t0\t%d\\n",\$1,\$2);}' '${fasta}.fai' |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter.bed


	if ${!exclude_bed.name.equals("NO_FILE")} ; then
		bedtools subtract -a TMP/jeter.bed -b "${exclude_bed}" |\
		awk -F '\t' 'int(\$2) < int(\$3)' |\
		sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter2.bed

		mv -v TMP/jeter2.bed TMP/jeter.bed
	fi

	mv -v TMP/jeter.bed select.bed
	test -s select.bed
	"""
	}

process SELECT01 {
	tag "${row.selname}"
	afterScript "rm -rf TMP"
	input:
		val(genomeId)
		val(row)
	output:
		tuple val(row),path("select.bed"),emit:output
	script:
		def genome = params.genomes[genomeId]
		def fasta = genome.fasta
		def gtf = genome.gtf
		def selid = row.selid
		def selname = row.selname
		def slop = row.slop
		def selexpr = row.selexpr
		def exclude_bed = row.exclude_bed==null?file("NO_FILE"):row.exclude_bed
	"""
	hostname 1>&2
	${moduleLoad("bedtools htslib")}
	set -o pipefail
	gunzip -c "${gtf}" |\\
		awk -F '\t' '${selexpr} {printf("%s\t%d\t%s\\n",\$1,int(\$4)-1,\$5);}' |\\
		bedtools slop -b ${slop} -g "${fasta}.fai" |\\
		sort -T . -t '\t' -k1,1 -k2,2n |\\
		bedtools merge > jeter.bed

	if ${!exclude_bed.name.equals("NO_FILE")} ; then

		bedtools subtract -a jeter.bed -b "${exclude_bed}" |\\
		awk -F '\t' 'int(\$2) < int(\$3)' |\\
		sort -T . -t '\t' -k1,1 -k2,2n > jeter2.bed

		mv -v jeter2.bed jeter.bed
	fi

	mv -v jeter.bed select.bed
	test -s select.bed
	"""
	}


process FIRST_INTRON {
	tag "${row.selid}"
	afterScript "rm -rf TMP"
	input:
		val(genomeId)
		val(row)
	output:
		tuple val(row),path("select.bed"),emit:output
	script:
		def genome = params.genomes[genomeId]
		def fasta = genome.fasta
		def gtf = genome.gtf
		def exclude_bed = row.exclude_bed==null?file("NO_FILE"):row.exclude_bed
	"""
	hostname 1>&2
	mkdir -p TMP
	${moduleLoad("bedtools")}
	
cat << __EOF__ > TMP/jeter.py
import sys
class Transcript:
    def __init__(self, chrom,strand, transcript_id):
        self.chrom = chrom
        self.strand = strand
        self.transcript_id = transcript_id
        self.exons = []

    def print_first_intron(self):
        # Trier les exons par position
        sorted_exons = sorted(self.exons, key=lambda exon: exon[0])

        # Si le transcript a au moins deux exons, calculer l'intron entre le premier et le deuxième exon
        if len(sorted_exons) >= 2:
            # En fonction de la direction (strand), retourner les coordonnées de l'intron
            if self.strand == '+':
                intron_start, intron_end = sorted_exons[0][1] +1 ,sorted_exons[1][0] - 1
            elif self.strand == '-':
                intron_start, intron_end = sorted_exons[-2][1]+1 ,sorted_exons[-1][0]- 1
            else:
                raise ValueError("Invalid value for 'strand'. Must be '+' or '-'.")

            print(f"{self.chrom}\t{intron_start -1}\t{intron_end}")

def filter_gtf_by_type():
    transcripts = {}

    for line in sys.stdin:
        if line.startswith('#'):
            continue

        fields = line.strip().split('\\t')
        if len(fields) < 9:
            continue

        feature = fields[2]
        if feature == "exon":
            transcript_id = None
            strand = fields[6]  # La colonne 7 (index 6) contient l'information sur le strand

            attributes = fields[8].split(';')
            for attr in attributes:
                if 'transcript_id' in attr:
                    transcript_id = attr.split('"')[1]

            if transcript_id:
                if transcript_id not in transcripts:
                    transcripts[transcript_id] = Transcript(fields[0],strand,transcript_id)

                start, end = int(fields[3]), int(fields[4])
                transcripts[transcript_id].exons.append((start, end))

    return transcripts

# Exemple d'utilisation
exon_transcripts = filter_gtf_by_type()

for tr   in exon_transcripts.values():
    tr.print_first_intron()
__EOF__


gunzip -c "${gtf}" | python3 |\
	sort -T TMP -t '\t' -k1,1 -k2,2n |\\
	bedtools merge > TMP/intron0.bed



	gunzip -c "${gtf}" |\\
		awk -F '\t' '(\$3=="exon") { printf("%s\t%d\t%s\\n",\$1,int(\$4)-1,\$5);}' |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\\
		bedtools merge > TMP/exons.bed

	bedtools subtract -a TMP/intron0.bed -b TMP/exons.bed |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\\
		awk -F '\t' 'int(\$2) < int(\$3)' > TMP/jeter.bed


	if ${!exclude_bed.name.equals("NO_FILE")} ; then

		bedtools subtract -a TMP/jeter.bed -b "${exclude_bed}" |\\
		awk -F '\t' 'int(\$2) < int(\$3)' |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter2.bed

		mv -v TMP/jeter2.bed TMP/jeter.bed
	fi

	mv -v TMP/jeter.bed select.bed
	test -s select.bed


	"""
	}

process INTRONS {
	tag "${row.selid}"
	afterScript "rm -rf TMP"
	input:
		val(genomeId)
		val(row)
	output:
		tuple val(row),path("select.bed"),emit:output
	script:
		def genome = params.genomes[genomeId]
		def fasta = genome.fasta
		def gtf = genome.gtf
		def exclude_bed = row.exclude_bed==null?file("NO_FILE"):row.exclude_bed
	"""
	hostname 1>&2
	${moduleLoad("bedtools")}
	set -o pipefail
	mkdir -p TMP

	gunzip -c "${gtf}" |\\
		awk -F '\t' '(\$3=="gene") { printf("%s\t%d\t%s\\n",\$1,int(\$4)-1,\$5);}' |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\\
		bedtools merge > TMP/genes.bed

	gunzip -c "${gtf}" |\\
		awk -F '\t' '(\$3=="exon") { printf("%s\t%d\t%s\\n",\$1,int(\$4)-1,\$5);}' |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\\
		bedtools merge > TMP/exons.bed

	bedtools subtract -a TMP/genes.bed -b TMP/exons.bed |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\\
		awk -F '\t' 'int(\$2) < int(\$3)' > TMP/jeter.bed


	if ${!exclude_bed.name.equals("NO_FILE")} ; then

		bedtools subtract -a TMP/jeter.bed -b "${exclude_bed}" |\\
		awk -F '\t' 'int(\$2) < int(\$3)' |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter2.bed

		mv -v TMP/jeter2.bed TMP/jeter.bed
	fi

	mv -v TMP/jeter.bed select.bed
	test -s select.bed
	"""
	}



process INTERGENIC {
	tag "${row.selid}"
	afterScript "rm -rf TMP"
	input:
		val(genomeId)
		val(row)
	output:
		tuple val(row),path("select.bed"),emit:output
	script:
		def genome = params.genomes[genomeId]
		def fasta = genome.fasta
		def gtf = genome.gtf
		def slop = row.slop
		def exclude_bed = row.exclude_bed==null?file("NO_FILE"):row.exclude_bed
	"""
	hostname 1>&2
	${moduleLoad("bedtools")}
	set -o pipefail
	mkdir -p TMP

	awk -F '\t' '{printf("%s\t0\t%d\\n",\$1,\$2);}' '${fasta}.fai' |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/genome.bed


	gunzip -c "${gtf}" |\\
		awk -F '\t' '(\$3=="gene") { printf("%s\t%d\t%s\\n",\$1,int(\$4)-1,\$5);}' |\\
		bedtools slop -b ${slop} -g "${fasta}.fai" |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\\
		bedtools merge > TMP/genes.bed

	bedtools subtract -a TMP/genome.bed -b TMP/genes.bed |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\\
		awk -F '\t' 'int(\$2) < int(\$3)' > TMP/jeter.bed


	if ${!exclude_bed.name.equals("NO_FILE")} ; then

		bedtools subtract -a TMP/jeter.bed -b "${exclude_bed}" |\\
		awk -F '\t' 'int(\$2) < int(\$3)' |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter2.bed

		mv -v TMP/jeter2.bed TMP/jeter.bed
	fi

	mv -v TMP/jeter.bed select.bed
	test -s select.bed
	"""
	}



process REGULATORY {
	tag "${row.selid}"
	afterScript "rm -rf TMP"
	input:
		val(bed)
		val(row)
	output:
		tuple val(row),path("select.bed"),emit:output
	script:
		def exclude_bed = row.exclude_bed==null?file("NO_FILE"):row.exclude_bed
	"""
	hostname 1>&2
	${moduleLoad("bedtools")}
	set -o pipefail
	mkdir -p TMP

	gunzip -c "${bed}" |\\
		awk -F '\t' '(\$4=="${row.selid}" || "${row.selid}" == "ALL_REG" )' |\\
		cut -f 1,2,3 |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\\
		bedtools merge > TMP/jeter.bed

	if ${!exclude_bed.name.equals("NO_FILE")} ; then

		bedtools subtract -a TMP/jeter.bed -b "${exclude_bed}" |\\
		awk -F '\t' 'int(\$2) < int(\$3)' |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter2.bed

		mv -v TMP/jeter2.bed TMP/jeter.bed
	fi

	mv -v TMP/jeter.bed select.bed
	test -s select.bed
	"""
	}


/** upstream / downstream gene */
process XXSTREAM {
	tag "${row.selid}"
	afterScript "rm -rf TMP"
	input:
		val(genomeId)
		val(row)
	output:
		tuple val(row),path("select.bed"),emit:output
	script:
		def genome = params.genomes[genomeId]
		def fasta = genome.fasta
		def gtf = genome.gtf
		def extend = row.extend
		def side = row.side
		def exclude_bed = (row.exclude_bed==null?file("NO_FILE"):row.exclude_bed)
		def left = "((C==\"-\" && ${side.equals("downstream")?1:0}) || (C==\"+\" && ${side.equals("upstream")?1:0}))"
	"""
	hostname 1>&2
	${moduleLoad("bedtools")}
	set -o pipefail
	mkdir -p TMP

	gunzip -c "${gtf}" |\\
		awk -F '\t' '(\$3=="gene") { printf("%s\t%d\t%s\t%s\\n",\$1,int(\$4)-1,\$5,\$7);}' |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/genes.strand.bed
	
	cut -f1,2,3 TMP/genes.strand.bed |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\\
		bedtools merge > TMP/genes.bed

	awk -F '\t' '{X=${extend};B=int(\$2);E=int(\$3);S=\$4; if(${left}) {E=E-1;B=E-X;if(B<0) B=0; if(E<0) E=0;} else {B=E;E=B+X;}  printf("%s\t%d\t%d\\n",\$1,B,E);}' TMP/genes.strand.bed |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\\
		awk -F '\t' 'int(\$2) < int(\$3)' |\\
		bedtools merge > TMP/xxx.stream.bed
		
	bedtools subtract -a TMP/xxx.stream.bed -b TMP/genes.bed |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\\
		awk -F '\t' 'int(\$2) < int(\$3)' > TMP/jeter.bed


	if ${!exclude_bed.name.equals("NO_FILE")} ; then

		bedtools subtract -a TMP/jeter.bed -b "${exclude_bed}" |\\
		awk -F '\t' 'int(\$2) < int(\$3)' |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter2.bed

		mv -v TMP/jeter2.bed TMP/jeter.bed
	fi

	mv -v TMP/jeter.bed select.bed
	test -s select.bed
	"""
	}


process DOWNLOAD_GFF_REG {
	input:
		val(genomeId)
	output:
		path("reg.gff.gz"),emit:output
		path("types.txt"),emit:types
	script:
		def genome = params.genomes[genomeId]
		def reference = genome.fasta
		def u = genome.ensembl_regulatory_gff_url
	"""
	hostname 1>&2
	${moduleLoad("jvarkit")}
	mkdir -p TMP
	set -o pipefail
	wget -O - "${u}" | gunzip -c |\\
		grep -v '^#' |\\
		java -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${reference}" --column 1 --convert SKIP |\\
		awk -F '\t' '{printf("%s\t%d\t%s\t%s\\n",\$1,int(\$4)-1,\$5,\$3);}' |\\
		LC_ALL=C sort -t '\t' -T TMP -k1,1 -k2,2n | uniq > reg.gff
	test -s reg.gff
	cut -f 4 reg.gff | sort | uniq > types.txt
	echo 'ALL_REG' >> types.txt
	test -s types.txt
	gzip reg.gff
	"""
}


process APPLY {
	tag "${row.selid} exclude:${row.exclude_bed}"
	afterScript "rm -rf TMP"
	errorStrategy "ignore"
	input:
		val(genomeId)
		path(vcf)
		path(cases)
		path(controls)
		val(row)
	output:
		path("${row.selid}.${params.step_id}.vcf.gz"),emit:vcf
		tuple val(row),path("${row.selid}.${params.step_id}.tsv"),emit:tsv
		tuple val(row),path("${row.selid}.${params.step_id}.exclude.bed"),emit:exclude
		path("${row.selid}.${params.step_id}.best_mqc.html"),emit:mqc
	when:
		row.selid.matches(params.selid_regex)	
	script:
		def genome = params.genomes[genomeId]
		def fasta = genome.fasta
		def gtf = genome.gtf
		def gnomad = genome.gnomad_genome
		def selid = row.selid
		def selbed = row.bed
	"""
	hostname 1>&2
	${moduleLoad("bedtools bcftools jvarkit")}
	set -o pipefail
	set -x
	mkdir -p TMP

	cat ${cases} ${controls} > TMP/all.samples.txt

	bcftools view --trim-alt-alleles --samples-file TMP/all.samples.txt  -O u --targets-file "${selbed}"  "${vcf}" |\\
	bcftools view -m2 -M 3 -O u   |\\
	bcftools norm -f "${fasta}" --multiallelics -any -O u |\\
	bcftools +contrast  -0 "${controls}" -1 "${cases}" -O v |\\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfgnomad --bufferSize 10000 --gnomad "${gnomad}" --fields "AF_POPMAX" |\\
	bcftools view  -i 'ALT!="*" && AC>0 && AF <= 0.05' -O z  -o TMP/mini.vcf.gz
	
	# TODO fix jvarkit
	touch TMP/best.properties
	
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar optimizefisher \\
		--threads ${task.cpus}  \\
		--duration '${params.duration}' \\
		--cases '${cases}' \\
		--controls '${controls}' -o TMP TMP/mini.vcf.gz 2> fisher.log
	
	# first and last line
	head -n1 TMP/best.tsv > TMP/best2.tsv
	tail -n 1 TMP/best.tsv >> TMP/best2.tsv

	touch TMP/exclude2.bed
	if ${row.exclude_bed!=null} ; then
		cat "${row.exclude_bed}" >> TMP/exclude2.bed
	fi

	# best hist as table

cat << __EOF__ > TMP/jeter.html
<!--
parent_id: ${params.step_id}
parent_name: "${params.step_name}"
parent_description: "Genetic Algorithm to find best set of parameters to optimize fisher test. VCF was <code>${vcf}</code>"
id: '${params.step_id}_${selid}_best'
section_name: '${row.selname}'
description: '${row.selname} ${params.step_name} : Best parameters.'
-->
__EOF__

	awk '(NR==1){for(i=0;i<=NF;i++) if(\$i=="INTERVAL") C=i;next;} {print \$C;}' TMP/best2.tsv |\\
		awk -F '[:-]' '{printf("<div>Best interval : <a target=\\"_blank\\" href=\\"https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=${genome.ucsc_name}&position=%s%%3A%s-%s\\">%s</a></div>\\n",\$1,int(\$2),\$3,\$0)}' >> TMP/jeter.html


	echo '<table class="table">' >> TMP/jeter.html
	paste	<(head -n 1 TMP/best2.tsv | tr "\t" "\\n") \\
		<(tail -n 1 TMP/best2.tsv | tr "\t" "\\n") |\\
		sed 's/&/\\&amp;/g; s/</\\&lt;/g; s/>/\\&gt;/g; s/"/\\&quot;/g; s/'"'"'/\\&#39;/g' |\\
		awk -F '\t' '{printf("<tr><th>%s</th><td>%s</td></tr>\\n",\$1,\$2);}'  >> TMP/jeter.html
	echo "</table>" >> TMP/jeter.html

	mv TMP/jeter.html "${selid}.${params.step_id}.best_mqc.html"
	
	# append best hit to exclude bed
	awk '(NR==1){for(i=0;i<=NF;i++) if(\$i=="INTERVAL") C=i;next;} {print \$C;}' TMP/best2.tsv |\\
		awk -F '[:-]' '{printf("%s\t%d\t%s\\n",\$1,int(\$2)-1,\$3)}' >> TMP/exclude2.bed

	
	mv TMP/best.vcf.gz ${selid}.${params.step_id}.vcf.gz
	mv TMP/best.tsv ${selid}.${params.step_id}.tsv
	mv TMP/best.properties ${selid}.${params.step_id}.properties

	sort -T TMP -t '\t' -k1,1 -k2,2n  "TMP/exclude2.bed" > "${selid}.${params.step_id}.exclude.bed"
	"""
	}


process MQC_TABLE {
tag "N=${L.size()}"
executor "local"
input:
	val(genomeId)
	val(L)
output:
	path("table_mqc.html"),emit:mqc
script:
	def genome = params.genomes[genomeId]
"""
mkdir -p TMP
export LC_ALL=C


cat << EOF > TMP/jeter.tsv
${L.join("\n")}
EOF


cat << EOF > TMP/jeter.txt
<!--
id: "table_for_${params.step_id}"
parent_id: "${params.step_id}"
parent_name: "${params.step_name}"
section_name: "Summary for ${params.step_name}"
description: "Summary for ${params.step_name}"
-->
<table class='table'><thead><caption>Summary ${params.step_name}. SOME COLUMNS MAY BE MISSING/SHIFTED BECAUSE JVARKIT REMOVES PARAMS IF THERE IS NO VARIANCE; NEED TO BE FIXED. </caption><tr><th>id</th>
EOF

# print header
head -n1 TMP/jeter.tsv | while IFS="\t" read N F
do
	head -n 1 "\${F}" |\
	tr "\t" "\\n" |\
	sed 's/&/\\&amp;/g; s/</\\&lt;/g; s/>/\\&gt;/g; s/"/\\&quot;/g; s/'"'"'/\\&#39;/g' |\\
	awk '{printf("<th>%s</th>\\n",\$0);}' >> TMP/jeter.txt
done
echo "</tr></thead><tbody>" >> TMP/jeter.txt

cat TMP/jeter.tsv | sort -t '\t' -k1,1V -T TMP| while IFS="\t" read N F
do
	# header and last line for best score
	head -n1 "\${F}" >  TMP/best2.tsv
	tail -n1 "\${F}" >> TMP/best2.tsv

	awk  -F '\t' -vN="\${N}" 'BEGIN{printf("<tr><td>%s</td>",N);} (NR==1){for(i=0;i<=NF;i++) if(\$i=="INTERVAL") C=i;next;} {for(i=1;i<=NF;i++) {S=\$i;if(i==C) { S=sprintf("<a target=\\"_blank\\" href=\\"https://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db=${genome.ucsc_name}&position=%s\\">%s</a>",S,S);} else {gsub(/^.*[<>]= /,"",S);} printf("<td>%s</td>",S);}} END {printf("</tr>\\n");}' TMP/best2.tsv  >> TMP/jeter.txt

done

echo "</tbody></table>" >> TMP/jeter.txt

mv TMP/jeter.txt table_mqc.html
"""
}
