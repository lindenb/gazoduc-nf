include {moduleLoad} from '../../modules/utils/functions.nf'

workflow GENETICALGO {
	take:
		genomeId
        	vcf
		cases
		controls
		exclude_bed_ch
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
				join(exclude_bed_ch.map{T->[T.selid,T]},remainder:true).
				map{T->T[1].plus("exclude_bed":T[2])}
			)

		intervals_ch = intervals_ch.mix(sel0b_ch.output)


		/** Whole genome */
		sel1a_ch = Channel.of([selid:"whole",selname:"All"])
	
		sel1b_ch  = WHOLE_GENOME(
			genomeId,
			sel1a_ch.map{T->[T.selid,T]}.
				join(exclude_bed_ch.map{T->[T.selid,T]},remainder:true).
				map{T->T[1].plus("exclude_bed":T[2])}
			)
		
		intervals_ch = intervals_ch.mix(sel1b_ch.output)

		/** introns */
		sel2a_ch = Channel.of([selid:"introns",selname:"Introns"])
		sel2b_ch  = INTRONS(
			genomeId,
			sel2a_ch.map{T->[T.selid,T]}.
				join(exclude_bed_ch.map{T->[T.selid,T]},remainder:true).
				map{T->T[1].plus("exclude_bed":T[2])}
			)

		intervals_ch = intervals_ch.mix(sel2b_ch.output)

		/** gene upstream/downstream */
		sel3a_ch = Channel.of(
			["upstream1kb","Upstream 1Kb","upstream",1_000],
			["upstream10kb","Upstream 10Kb","upstream",10_000],
			["upstream100kb","Upstream 100Kb","upstream",100_000],
			["upstream1Mb","Upstream 1M","upstream",1_000_000],
                        ["downstream1kb","Downstream 1Kb","downstream",1_000],
                        ["downstream10kb","Downstream 10Kb","downstream",10_000],
                        ["downstream100kb","Downstream 100Kb","downstream",100_000],
                        ["downstream1Mb","Downstream 1M","downstream",1_000_000]
			).
			map{T->[ selid:T[0], selname:T[1], side:T[2], extend:T[3]]}

		sel3b_ch = XXSTREAM (
			genomeId,
			sel3a_ch.map{T->[T.selid,T]}.
				join(exclude_bed_ch.map{T->[T.selid,T]},remainder:true).
				map{T->T[1].plus("exclude_bed":T[2])}
			)
		intervals_ch = intervals_ch.mix(sel3b_ch.output)

		APPLY(  genomeId,
			vcf,
			cases,
			controls,
			intervals_ch.map{T->T[0].plus(bed:T[1])}
			)
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
	awk -F '\t' '{printf("%s\t0\t%d\\n",\$1,\$2);}' '${fasta}.fai' > TMP/jeter.bed


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
		awk -F '\t' '(\$3=="exons") { printf("%s\t%d\t%s\\n",\$1,int(\$4)-1,\$5);}' |\\
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


process APPLY {
	tag "${row.selid}"
	afterScript "rm -rf TMP"
	cpus 10
	memory "10g"
	input:
		val(genomeId)
		path(vcf)
		path(cases)
		path(controls)
		val(row)
	output:
		path("${row.selid}.${params.step_id}.vcf.gz"),emit:vcf
		path("${row.selid}.${params.step_id}.tsv"),emit:tsv
		tuple val(row),path("${selid}.${params.step_id}.exclude.bed"),emit:exclude
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
		--cases '${cases}' \\
		--controls '${controls}' -o TMP TMP/mini.vcf.gz
	
	# first and last line
	head -n1 TMP/best.tsv > TMP/best2.tsv
	tail -n 1 TMP/best.tsv >> TMP/best2.tsv

	touch TMP/exclude2.bed
	if ${row.exclude_bed==null} ; then
		cat "${row.exclude_bed}" >> TMP/exclude2.bed
	fi

	# best hist as table
	echo "<pre>" > TMP/jeter.html
	paste	<(head -n 1 | TMP/best2.tsv | tr "\t" "\\n") \\
		<(tail -n 1 | TMP/best2.tsv | tr "\t" "\\n") |\\
		column -t -s '\t'  >> TMP/jeter.html
	echo "</pre>" > TMP/jeter.html

	# append best hit to exclude bed
	awk '(NR==1){for(i=0;i<=NF;i++) if(\$i=="INTERVAL") C=i;next;} {print \$C;}' TMP/best2.tsv |\\
		awk -F '[:-]' '{printf("%s\t%d\t%s\\n",\$1,int(\$2)-1,\$3)}' >> TMP/exclude2.bed
	
	mv TMP/best.vcf.gz ${selid}.${params.step_id}.vcf.gz
	mv TMP/best.tsv ${selid}.${params.step_id}.tsv
	mv TMP/best.properties ${selid}.${params.step_id}.properties

	mv "TMP/exclude2.bed" "${selid}.${params.step_id}.exclude.bed"
	"""
	}


