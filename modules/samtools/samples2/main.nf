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
process SAMTOOLS_SAMPLES {
tag "${meta.id?:""}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path("REFS/*")
    tuple val(meta ),path(samplesheet)
output:
    tuple val(meta),path("*.csv"),emit:samplesheet
    path("versions.yml"),emit:versions
script:
    def args1 = task.ext.args1?:""
    def prefix = task.ext.prefix?"${meta.id}"
"""
mkdir -p TMP

cat << 'EOF' | R --vanilla
data <- read.csv("${samplesheet}", stringsAsFactors = FALSE)

# Check if required columns exist
if (!("bam" %in% names(data))) {
  stop("Missing required column 'bam' in ${samplesheet}.")
}
filtered_data <- data[data\$bam != "" & data\$bam != ".", ]
# Write the 'bam' column to a text file
writeLines(filtered_data\$bam, "TMP/bams.list")
EOF

cat  TMP/bams.list |\\
    samtools samples  \\
        `find REFS/ \\( -name "*.fasta" -o  -name "*.fa" -o -name "*.fna" \\)  -printf "-f %p "`  > TMP/samplesheet2.tsv


cat << 'EOF' | R --vanilla
data1 <- read.csv("${samplesheet}", stringsAsFactors = FALSE)

# Step 2: Check required columns in data1
if (!all(c("sample", "bam") %in% names(data1))) {
  stop("Missing required columns 'sample' or 'bam' in data1")
}

# Step 3: Read data2 (tab-delimited, no header)
data2_raw <- read.table("data2.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(data2_raw) <- c("sample", "bam", "fasta")

# Step 4: Merge data1 and data2 on 'bam'
merged <- merge(data1, data2_raw, by = "bam", suffixes = c(".data1", ".data2"), all.x = TRUE)

# Step 5: Replace empty or '.' in data1$sample
replace_sample <- merged$sample.data1 == "" | merged$sample.data1 == "."
merged$sample.data1[replace_sample] <- merged$sample.data2[replace_sample]

# Step 6: Handle 'fasta' column if it exists in data1
if ("fasta" %in% names(data1)) {
  # Replace empty or '.' fasta in data1 with data2's fasta
  replace_fasta <- merged$fasta == "" | merged$fasta == "."
  merged$fasta[replace_fasta] <- merged$fasta.data2[replace_fasta]

  # Check for conflicts: both non-empty and not equal
  conflict <- !(merged$fasta == "" | merged$fasta == ".") &
              !(merged$fasta.data2 == "" | merged$fasta.data2 == ".") &
              merged$fasta != merged$fasta.data2

  if (any(conflict, na.rm = TRUE)) {
    stop("Conflict detected: fasta values differ between data1 and data2 for some rows.")
  }

} else {
  # If 'fasta' not in data1, just use fasta from data2
  merged$fasta <- merged$fasta.data2
}

# Step 7: Construct final output data frame
# Replace original sample with updated one
merged$sample <- merged$sample.data1

# Drop intermediate columns
final <- merged[, c("sample", "bam", "fasta")]

# Step 8: Write final merged data to CSV with header
write.csv(final, "merged_output.csv", row.names = FALSE)
EOF


mv TMP/samplesheet2.tsv ./

cat << END_VERSIONS > versions.yml
"${task.process}":
	samtools: "\$(samtools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""

stub:
"""
find \${PWD}/BAMS/ -name "*am" | sort -V | awk '{printf("S%d\t%s\t%s\\n",NR,\$0,"${params.fasta?:"."}");}'  > samplesheet.tsv
touch versions.yml
"""
}
