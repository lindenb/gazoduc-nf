/**

original workflow was designed by  Molitor Corentin  https://github.com/MCorentin/kilda


                    GNU AFFERO GENERAL PUBLIC LICENSE
                       Version 3, 19 November 2007

 Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
 Everyone is permitted to copy and distribute verbatim copies
 of this license document, but changing it is not allowed.

                            Preamble

  The GNU Affero General Public License is a free, copyleft license for
software and other kinds of works, specifically designed to ensure
cooperation with the community in the case of network server software.

  The licenses for most software and other practical works are designed
to take away your freedom to share and change the works.  By contrast,
our General Public Licenses are intended to guarantee your freedom to
share and change all versions of a program--to make sure it remains free
software for all its users.

  When we speak of free software, we are referring to freedom, not
price.  Our General Public Licenses are designed to make sure that you
have the freedom to distribute copies of free software (and charge for
them if you wish), that you receive source code or can get it if you
want it, that you can change the software or use pieces of it in new
free programs, and that you know you can do these things.

  Developers that use our General Public Licenses protect your rights
with two steps: (1) assert copyright on the software, and (2) offer
you this License which gives you legal permission to copy, distribute
and/or modify the software.

  A secondary benefit of defending all users' freedom is that
improvements made in alternate versions of the program, if they
receive widespread use, become available for other developers to
incorporate.  Many developers of free software are heartened and
encouraged by the resulting cooperation.  However, in the case of
software used on network servers, this result may fail to come about.
The GNU General Public License permits making a modified version and
letting the public access it on a server without ever releasing its
source code to the public.

  The GNU Affero General Public License is designed specifically to
ensure that, in such cases, the modified source code becomes available
to the community.  It requires the operator of a network server to
provide the source code of the modified version running there to the
users of that server.  Therefore, public use of a modified version, on
a publicly accessible server, gives the public access to the source
code of the modified version.

  An older license, called the Affero General Public License and
published by Affero, was designed to accomplish similar goals.  This is
a different license, not a version of the Affero GPL, but Affero has
released a new version of the Affero GPL which permits relicensing under
this license.

*/


workflow {
	def fasta = file(params.fasta)
	def fai = file(params.fasta+".fai")

	data_ch = PREPARE_DATA(fasta,fai)

	ch1 =  data_ch.KIV2_bed.map{["kiv2",it,6]}.mix(data_ch.LPA_bed.map{["norm",it,1]})

	ch2 = PREPARE_KMERS(fasta, fai, ch1)

	ch3 = ch2.output.combine(ch2.output).
		filter{it[0].equals("kiv2") && it[2].equals("norm")}.
		map{[it[1],it[3]]}

	ch4 = REMOVE_COMMON_KMERS(ch3)

	ch5 = OUTPUT_FASTA(ch4.kiv2.combine(ch4.norm))

	ch6 = CREATE_FASTA_KMERS(ch4.kiv2.combine(ch4.norm).combine(data_ch.lpa_three_rsids_tsv))

	ch7 = Channel.fromPath(params.samplesheet).splitCsv(header:true,sep:'\t').
		map{[it.sample,it.bam,it.bai]}.
		combine(data_ch.merged_bed).combine(ch6.output)

	ch8 = COUNT_KMERS(fasta,fai,ch7)

	ch9 = GET_SCRIPT()

	ch10 = KILDA(ch4.kiv2,ch4.norm,data_ch.lpa_three_rsids_tsv,ch9.output, ch8.output.collect())

	ZIPIT(ch10.output)
	}

process PREPARE_DATA {
input:
	path(fasta)
	path(fai)
output:
	path("KIV2.bed"),emit:KIV2_bed
	path("LPA.bed"),emit:LPA_bed
	path("merged.bed"),emit:merged_bed
	path("lpa_rsids.tsv"),emit:lpa_rsids_tsv
	path("lpa_three_rsids.tsv"),emit:lpa_three_rsids_tsv
script:
	def build= (fasta.name.contains("hs37d5")?"hg19":(fasta.name.contains("hs38me")?"hg38":params.build))

"""

if ${build.equals("hg38")} ; then

echo -e 'chr6\t160611534\t160646545' > KIV2.bed

echo -e 'chr6\t160531483\t160611532\\nchr6\t160646546\t160666375' > LPA.bed

elif ${build.equals("hg19")} ; then


echo -e 'chr6\t161032566\t161067577' > KIV2.bed

echo -e "chr6\t160952515\t161032564\\nchr6\t161067578\t161087407" > LPA.bed

else
	echo "BOUM ${build}" 1>&2
fi

cat KIV2.bed LPA.bed | sort -t '\t' -k1,1 -k2,2n | bedtools merge > merged.bed

cat << EOF > lpa_rsids.tsv
rs41272114	AACATGATAGACATACGCATTTGGATAGTAT	AACATGATAGACATATGCATTTGGATAGTAT
rs10455872	TGTTCTCAGAACCCAATGTGTTTATACAGGT	TGTTCTCAGAACCCAGTGTGTTTATACAGGT
rs3798220	CAGCCTAGACACTTCTATTTCCTGAACATGA	CAGCCTAGACACTTCCATTTCCTGAACATGA
rs186696265	GACCACGACGAACGACGACACGAGCAAACCA	GACCACGACGAACGATGACACGAGCAAACCA
rs1623955	AGTGCCCAGAAAGTGTGTCCCAATCCCAGGA	AGTGCCCAGAAAGTGGGTCCCAATCCCAGGA
rs75692336	TGCTACATAAGCCTGCTGATTCCAAAGTCTT	TGCTACATAAGCCTGATGATTCCAAAGTCTT
rs41259144	ACAGGATCTGGATTTCGGCAGTAGTTCTTGA	ACAGGATCTGGATTTTGGCAGTAGTTCTTGA
rs140570886	AATTTAACTATGATGTACCCCAGTAAATCAG	AATTTAACTATGATGCACCCCAGTAAATCAG
rs41267813	CTGGGGTCCTCTGATGCCAGTGTGGTATCAT	CTGGGGTCCTCTGATACCAGTGTGGTATCAT
rs41267811	GTCCCTTCTGTGTCTGAGCATCGCGTCAGGT	GTCCCTTCTGTGTCTCAGCATCGCGTCAGGT
rs139145675	CACCATCAGGGTTACGGCAGTACTGAAAACA	CACCATCAGGGTTACAGCAGTACTGAAAACA
rs41267809	TCATTCTCAATAACAAGGAGCTGGGCTTCCT	TCATTCTCAATAACAGGGAGCTGGGCTTCCT
rs3124784	CAGGCTTATTGGGGCGTGCACAGCCAAGACC	CAGGCTTATTGGGGCATGCACAGCCAAGACC
rs7412	GAAATATTTGGTTTACAGAAGTACTCACTGG	GAAATATTTGGTTTATAGAAGTACTCACTGG
EOF



cat << EOF > lpa_three_rsids.tsv
rs41272114	AACATGATAGACATACGCATTTGGATAGTAT	AACATGATAGACATATGCATTTGGATAGTAT
rs10455872	TGTTCTCAGAACCCAATGTGTTTATACAGGT	TGTTCTCAGAACCCAGTGTGTTTATACAGGT
rs3798220	CAGCCTAGACACTTCTATTTCCTGAACATGA	CAGCCTAGACACTTCCATTTCCTGAACATGA
EOF


"""
}

workflow PREPARE_KMERS {
take:
	fasta
	fai
	row
main:
	ch1 = row.map{it.plus(true)}.mix(row.map{it.plus(false)})
	ch2 = COUNT_KMERS_REGION(fasta ,fai ,ch1 )
	ch3 = FILTER_ON_OCCURENCE(ch2.output.filter{it[1]==false})


	ch3 = ch3.output.combine(ch2.output.filter{it[1]==true}).
		filter{it[0].equals(it[2])}.
		map{[it[0],it[1],it[5]]}



	ch5 = FILTER_KMERS_OCCURING_OUTSIDE_REGION(ch3)

emit:
	output = ch5.output

}




process COUNT_KMERS_REGION {
    tag "${type} ${inverse}"
    label "process_medium"
    afterScript "rm -rf TMP"
    input:
        path(fasta)
        path(fai)
	tuple val(type),path(region_bed),val(occ),val(inverse)
    output:
    	    tuple val(type),val(inverse),val(occ),path("${type}${inverse?".complement":""}.kmers"),emit:output
    script:        
        """
	set -o pipefail
	mkdir -p TMP
	set -x

	if ${inverse} ; then

	        awk -F '\t' '{ print \$1"\\t0\\t"\$2 }' '${fai}' > TMP/genome.bed
        
        	bedtools intersect -v -a TMP/genome.bed -b ${region_bed} |\\
			awk -F '\t' '{printf("%s:%d-%s\\n",\$1,int(\$2)+1,\$3);}' > TMP/jeter.regions
	else
	
	        awk -F '\t' '{printf("%s:%d-%s\\n",\$1,int(\$2)+1,\$3);}' "${region_bed}" > TMP/jeter.regions
	
	fi

        samtools faidx "${fasta}" -r TMP/jeter.regions > TMP/jeter.fa

        jellyfish count -C -t ${task.cpus} -m ${params.kmer_size} -s 10G -o TMP/kmers.data TMP/jeter.fa

	mv TMP/kmers.data "${type}${inverse?".complement":""}.kmers"
	"""
}

process FILTER_ON_OCCURENCE {
    tag "${type}"
    label "process_medium"
    afterScript "rm -rf TMP"
    input:
	tuple val(type),val(inverse),val(occ),path(kmers)
    output:
	tuple val(type),path("${type}.${occ}copies.dump"),emit:output
    script:
	"""
        jellyfish dump -c -t -L ${occ} -U ${occ} -o "${type}.${occ}copies.dump" "${kmers}"
	"""
}

 
process FILTER_KMERS_OCCURING_OUTSIDE_REGION {
    tag "${type}"
    label "process_medium"
    afterScript "rm -rf TMP"
    input:
    	tuple val(type),  path(kmers_filtered_dump), path(kmers_outside_region_count)
    output:
        tuple val(type),path("${type}.specific.fasta"),emit:output
    script:
        kmers_filtered_fasta = "${kmers_filtered_dump.baseName}.fasta"
        kmers_specific = "${kmers_filtered_dump.baseName}_specific.fasta"
        
        """
        set -eo pipefail
        mkdir -p TMP
	set -x

        cut -f 1 '${kmers_filtered_dump}' | awk '{print ">kmer_"NR"\\n"\$1 }' > "TMP/jeter.fasta"

        jellyfish query -s "TMP/jeter.fasta" -o TMP/jeter.counts ${kmers_outside_region_count}
        
        awk '\$2 == 0 { print \$1 }' TMP/jeter.counts > TMP/specific_kmers.list

        cut -f 1 TMP/specific_kmers.list | LC_ALL=C sort -T TMP > "${type}.specific.fasta"
        """
}

process REMOVE_COMMON_KMERS {
    label "process_short"
    input:
        tuple path(kiv2_kmers_list),path(norm_kmers_list)        
    output:
        path("kiv2.tsv"),emit:kiv2
        path("norm.tsv"),emit:norm
    script:
    
        """
	LC_ALL=C comm -13 ${kiv2_kmers_list} ${norm_kmers_list} > norm.tsv
	LC_ALL=C comm -23 ${kiv2_kmers_list} ${norm_kmers_list} > kiv2.tsv
        """
}


process OUTPUT_FASTA {
    label "process_short"
    input:
        tuple(path(kiv2_unique_kmers_list), path(norm_unique_kmers_list))
    
    output:
        tuple(path(kiv2_unique_kmers_fasta), path(norm_unique_kmers_fasta))

        
    script:
        kiv2_unique_kmers_fasta = "${kiv2_unique_kmers_list.baseName}.fasta"
        norm_unique_kmers_fasta = "${norm_unique_kmers_list.baseName}.fasta"
        
        """
        set -eo pipefail
        
        awk '{print ">kmer_KIV2_"NR"\\n"\$1 }' ${kiv2_unique_kmers_list} > ${kiv2_unique_kmers_fasta}
        awk '{print ">kmer_NORM_"NR"\\n"\$1 }' ${norm_unique_kmers_list} > ${norm_unique_kmers_fasta}
        """
}


process CREATE_FASTA_KMERS {
    label "process_short"
    input:
       tuple path(kiv2_kmers), path(norm_kmers), path(rsids_list)
    
    output:
        path("NORM_KIV2_kmers.fasta"),emit:output
    script:
        """
        awk '{ print ">KIV2_"NR"\\n"\$1 }' ${kiv2_kmers} > jeter.fa
        awk '{ print ">NORM_"NR"\\n"\$1 }' ${norm_kmers} >> jeter.fa
        awk '{ print ">"\$1"_ref\\n"\$2"\\n>"\$1"_alt\\n"\$3 }' ${rsids_list} >> jeter.fa
	mv jeter.fa NORM_KIV2_kmers.fasta
        """
}

process COUNT_KMERS {
tag "${sample}"
label "process_short"
afterScript "rm -rf TMP"
input:
	path(fasta)
	path(fai)
	tuple val(sample),path(bam),path(bai),path(bed),path(kmers)
output:
	path("${sample}.counts"),emit:output
script:
"""
mkdir -p TMP
samtools view -F 3844  --uncompressed -O BAM  -M -L ${bed} --threads 1 --reference "${fasta}" "${bam}" |\\
samtools fastq > TMP/jeter.fq

jellyfish count -t ${task.cpus} -m ${params.kmer_size} -s 100M -C -o TMP/jeter.kmers --if=${kmers} TMP/jeter.fq
jellyfish dump TMP/jeter.kmers -c -t > TMP/jeter.counts

mv TMP/jeter.counts "${sample}.counts"
"""
}

process GET_SCRIPT {
executor "local"
output:
	path("kilda.py"),emit:output
script:
"""
wget -O kilda.py "https://raw.githubusercontent.com/MCorentin/kilda/refs/heads/main/bin/kilda.py"
"""
}

process KILDA {
label "process_medium"
input:
	path(kiv2_kmers)
        path(norm_kmers)
	path(rsid)
	path(kilda_py)
	path(counts)
output:
	path("OUTPUT"),emit:output
script:
"""
mkdir -p TMP

find . -name "*.counts" | sed 's/^\\.\\///' | awk -F '.' '{printf("%s\t%s\\n",\$1,\$0);}' > TMP/counts.list


python3 ${kilda_py} \
        -c TMP/counts.list -o OUTPUT -v -p \
        -k ${kiv2_kmers} -l ${norm_kmers} \
        ${rsid}

"""
}

process ZIPIT {
input:
	path(dir)
output:
	path("${params.prefix}.output.zip"),emit:output
script:
"""
zip -r "${params.prefix}.output.zip" "${dir}"
"""
}
