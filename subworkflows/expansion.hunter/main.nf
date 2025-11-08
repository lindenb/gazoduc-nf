/*

Copyright (c) 2025 Pierre Lindenbaum

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
include { XHUNTER_APPLY                } from '../../modules/expansion.hunter/apply'
include { BCFTOOLS_MERGE               } from '../../modules/bcftools/merge'
include { isBlank                      } from '../../modules/utils/functions.nf'
include {GHOSTSCRIPT_MERGE             } from '../../modules/gs/merge'

workflow EXPANSION_HUNTER {
	take:
		workflow_metadata
		fasta
        fai
        dict
        catalog
        bams
	main:
        versions = Channel.empty()

        XHUNTER_APPLY(
            fasta,
            fai,
            catalog,
            bams
            )
		versions = versions.mix(XHUNTER_APPLY.out.versions)

        BCFTOOLS_MERGE(
	    [[id:"nobed"],[]],
            XHUNTER_APPLY.out.vcf
                .map{meta,vcf,tbi->[vcf,tbi]}
                .flatMap()
                .collect()
                .map{[
                    [id:workflow_metadata.id],
                    it.sort()
                    ]}
            )
        versions = versions.mix(BCFTOOLS_MERGE.out.versions)
        
        sn2status_ch = bams.map{meta,_bam,_bai->meta}
            .filter{meta->!isBlank(meta.status) && meta.status.matches("(case|control)")}
            .map{meta->meta.id+"\t"+meta.status}
            .collectFile(name: 'sample2status.tsv', newLine: true)
            .map{fn->[[id:"xhunter"],fn]}
	
	what = Channel.of("average","max","min")

        FISHER_TEST(
            BCFTOOLS_MERGE.out.vcf.map{meta,vcf,tbi->[meta,vcf,tbi]},
            sn2status_ch
		.filter{_meta,f->(f instanceof java.nio.file.Path)}
		.combine(what)
		.map{meta,f,w->[meta.plus([what:w]),f]}
            )

	GHOSTSCRIPT_MERGE(
		FISHER_TEST.out.pdf
			.flatMap{meta,pdfs->{
				def L = (pdfs instanceof List?pdfs:[pdfs]);
				def L2 = [];
				for(int i=0;i< L.size();i++) {
					L2.add([[id: meta.what],L[i]]);
					}
				return L2;
				}}
			.groupTuple()
			.view()
		)
	versions = versions.mix(GHOSTSCRIPT_MERGE.out.versions)
	emit:
		versions
        vcf = BCFTOOLS_MERGE.out.vcf
	}

process FISHER_TEST {
tag "${meta.id} ${meta2.what}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.02.yml"
afterScript "rm -rf TMP"
input:
    tuple val(meta) ,path(vcf),path(tbi)
    tuple val(meta2),path(sample2status)
output:
    tuple val(meta2),path("*.genes.pdf"),emit:pdf
    tuple val(meta2),path("*.fisher.tsv"),emit:fisher
    path("versions.yml"),emit:versions
script:
    def prefix=task.ext.prefix?:"${meta.id}"
    def what = meta2.what
"""
mkdir -p TMP/PDFS TMP/BEDS

cat << 'EOF' > TMP/max.awk
	{
	j=0;
	M = 0.0;
	N = split(\$3,a,/[\\/]/);
	for(i=1;i<=N;i++) {
		if(a[i]=="./." || a[i]==".") continue;
		V= int(a[i])*1.0;
		if(j==0 || V < M) M=V;
		j++;
		}
	if(j==0) next;
	printf("%s\t%s\t%f\t%s\\n",\$1,\$2,M,\$4);
	}
EOF

cat << 'EOF' > TMP/min.awk
	{
	j=0;
	M = 0.0;
	N = split(\$3,a,/[\\/]/);
	for(i=1;i<=N;i++) {
		j=1;
		if(a[i]=="./." || a[i]==".") continue;
		V= int(a[i])*1.0;
		if(j==0 || V > M) M=V;
		j++;
		}
	if(j==0) next;
	printf("%s\t%s\t%f\t%s\\n",\$1,\$2,M,\$4);
	}
EOF


cat << 'EOF' > TMP/average.awk
	{
	T = 0.0;
	N = split(\$3,a,/[\\/]/);
	j = 0;
	for(i=1;i<=N;i++) {
		if(a[i]=="./." || a[i]==".") continue;	
		V=int(a[i]);
		if(V==0) continue;
		T+=V;
		j++;
		}
	if(j==0) next;
	printf("%s\t%s\t%f\t%s\\n",\$1,\$2,T/j,\$4);
	}
EOF



i=1
bcftools query -f '%CHROM\t%POS0\t%END\\n' '${vcf}' |\\
	 split -a 9 --additional-suffix=.bed --lines=1000 - TMP/chunck.

find TMP/ -type f -name "chunck.*.bed" | while read BED
do

echo "\${BED}" 1>&2

bcftools view -O u --regions-file "\${BED}"  '${vcf}'  |\\
    bcftools query -f '[%SAMPLE\t%ID\t%REPCN\t%CHROM:%POS-%END\\n]' -  |\\
    awk -F '\t' -f TMP/${what}.awk |\\
    sort  -S ${task.memory.kilo} -t '\t' -T TMP -k1,1 > TMP/jeter.a

sort -t '\t'  -S ${task.memory.kilo} -T TMP -k1,1 '${sample2status}'  > TMP/jeter.b

join -t '\t' -1 1 -2 1 -o '1.1,1.2,1.3,1.4,2.2' TMP/jeter.a TMP/jeter.b > TMP/jeter.tsv


cat << '__EOF__' > TMP/jeter.R

dat <- read.table("TMP/jeter.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")

colnames(dat) <- c("sample", "gene", "average", "position","status")

# Convert average to numeric, status to factor
dat\$average <- as.numeric(dat\$average)

# Ensure status is factor 
 dat\$status <- factor(dat\$status, levels = c("case", "control"))

# Get unique genes and sort
genes <- sort(unique(dat\$gene))

# Open PDF and create one page per gene
pdf("TMP/genes.pdf", width = 6, height = 6)
for (g in genes) {
  d <- dat[dat\$gene == g, , drop = FALSE]
  if (nrow(d) == 0) next

  # require two status 
  if (length(unique(d\$status)) == 2) {
   boxplot(
        average ~ status,
        data = d,
        main = paste("Locus: ", g),
	sub = "${what}",
        ylab = "average",
        col = c("#66c2a5", "#fc8d62")
        ) # two distinct colors
    # overlay points for raw data
    stripchart(
        average ~ status, 
        data = d, 
        vertical = TRUE,
        add = TRUE,
        method = "jitter",
        pch = 21,
        col = "black",
        bg = "white"
        )
  }
}
dev.off()

res_df <- data.frame(
        position = character(),
        gene = character(),
        median_control = double(),
        above_control = integer(),
        below_control = integer(),
        above_case = integer(),
        below_case = integer(),
        odds_ratio = integer(),
        p_value = double(),
        stringsAsFactors = FALSE
        )

for (g in genes) {
    d <- dat[dat\$gene == g, , drop = FALSE]
    n_control <- sum(d\$status == "control", na.rm = TRUE)
    if (n_control == 0) next
    position_str <-  sort(unique(dat\$position[dat\$gene == g]))[1]


    median_control <- median(d\$average[d\$status == "control"], na.rm = TRUE)

    # classify >= median as TRUE; treat NA averages as NA (we'll not count them in cells)
    is_above <- ifelse(is.na(d\$average), NA, d\$average >= median_control)

    above_control <- sum(is_above & d\$status == "control", na.rm = TRUE)
    below_control <- sum((!is_above) & d\$status == "control", na.rm = TRUE)
    above_case <- sum(is_above & d\$status == "case", na.rm = TRUE)
    below_case <- sum((!is_above) & d\$status == "case", na.rm = TRUE)


    # Build contingency table: rows = >=m (above), <m (below); cols = control, case
    tab <- matrix(c(above_control, above_case, below_control, below_case),
        nrow = 2, byrow = TRUE,
        dimnames = list(c("ge_median", "lt_median"), c("control", "case")))
    
    ft <- fisher.test(tab)

    odds_ratio <- NA
    p_value <- ft\$p.value
    if (!is.null(ft\$estimate)) odds_ratio <- as.numeric(ft\$estimate)

    res_df[nrow(res_df) + 1,] = c(
        position_str,
        g,
        median_control,
        above_control,
        below_control,
        above_case,
        below_case,
        odds_ratio,
        p_value
        )
  }

write.table(res_df, file = "TMP/fisher.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
__EOF__


R --vanilla --quiet < TMP/jeter.R

if test -s TMP/genes.pdf
then
  mv TMP/genes.pdf "TMP/PDFS/\${i}.genes.pdf"
fi

if test -s TMP/fisher.tsv
then
  mv TMP/fisher.tsv "TMP/\${i}.fisher.tsv"
fi

i=\$((i+1))

done

mv TMP/PDFS/*.pdf ./

find ./TMP -name "*.fisher.tsv" -exec cat '{}' ';' |\\
	sed 's/^position/#position/' |\\
	LC_ALL=C sort -T TMP -t '\t' -k9,9g | uniq >  "${prefix}.${what}.fisher.tsv"

touch versions.yml
"""
}




