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
process PLOT_ASSOC {
label "process_single"
tag "${assoc.name}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta2),path(fai)
	tuple val(meta ),path(assoc)
output:
    tuple val(meta ),path("*.png"),optional:true,emit:png
    tuple val(meta ),path("*.pdf"),optional:true,emit:pdf
    path("versions.yml"),emit:versions
script:
    def format = task.ext.format?:"png"
    def prefix = task.ext.prefix?:"${assoc.baseName}" //not meta.id please
    def plot_size = task.ext.plot_size?:500
"""
hostname 1>&2
mkdir -p TMP
set -x

cat << '__EOF__' | R --no-save 
## R-base (no packages) Manhattan-like plot from PLINK output + .fai
## - Y axis is -log10(P)
## - Alternating point colors by chromosome
## - Background vertical bars (alternating) to delimit chromosomes

## Inputs

## Read PLINK output
d <- read.table("${assoc}", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
head(d)
d\$CHR <- as.character(d\$CHR)
d\$BP  <- as.integer(d\$BP)
d\$P   <- as.numeric(d\$P)

## Read FAI (only first two columns are needed: chr name, chr length)
fai <- read.table("${fai}", header=FALSE, sep="\t", stringsAsFactors=FALSE, quote="")
fai <- fai[, 1:2]
colnames(fai) <- c("CHR", "LEN")
fai\$CHR <- as.character(fai\$CHR)
fai\$LEN <- as.numeric(fai\$LEN)
head(fai)

## Normalize chromosome names: drop optional leading "chr"/"Chr"/"CHR"
fai\$CHR <- sub("^(chr)", "", fai\$CHR, ignore.case=TRUE)
d\$CHR   <- sub("^(chr)", "", d\$CHR,   ignore.case=TRUE)


## Keep only variants on chromosomes present in FAI (and in the same order as FAI)
d <- d[d\$CHR %in% fai\$CHR, , drop=FALSE]
d\$CHR <- factor(d\$CHR, levels=fai\$CHR, ordered=TRUE)
d <- d[order(d\$CHR, d\$BP), , drop=FALSE]





## Genome-wide coordinate: cumulative chromosome starts
chr_start <- c(0, cumsum(fai\$LEN))[seq_len(nrow(fai))]  # start position for each chr
names(chr_start) <- fai\$CHR

d\$X <- d\$BP + chr_start[as.character(d\$CHR)]

## Y: -log10(P). Protect against P==0.
d\$P[d\$P <= 0] <- min(d\$P[d\$P > 0], na.rm=TRUE) / 10
d\$Y <- -log10(d\$P)

## Axis ticks at chromosome midpoints
chr_mid <- chr_start + fai\$LEN / 2
names(chr_mid) <- fai\$CHR

## Colors
point_cols <- c("#1f77b4", "#ff7f0e")     # alternating point colors
bg_cols    <- c("grey95", "white")        # alternating chromosome background bands

## Plot

${format}("TMP/${prefix}_mqc.${format}",width = ${((plot_size as double) * 5.0 ) as int}, height = ${plot_size}, unit = "px")


op <- par(mar=c(5, 4, 3, 1) + 0.1)
on.exit(par(op), add=TRUE)

plot(d\$X, d\$Y,
     type="n",
     xlab="Genomic position (by chromosome)",
     ylab="-log[10](P)",
     xaxt="n",
     main="${assoc.name}"
     )

## Background chromosome bands + boundary lines
ymax <- max(d\$Y, na.rm=TRUE)
for (i in seq_len(nrow(fai))) {
  x0 <- chr_start[i]
  x1 <- chr_start[i] + fai\$LEN[i]
  rect(x0, 0, x1, ymax, col=bg_cols[(i %% 2) + 1], border=NA)
  abline(v=x0, col="grey80", lwd=1)  # boundary at chromosome start
}
abline(v=sum(fai\$LEN), col="grey80", lwd=1) # last boundary

## Points (alternating by chromosome)
chr_index <- as.integer(d\$CHR)
cols <- point_cols[(chr_index %% 2) + 1]
points(d\$X, d\$Y, pch=20, cex=0.6, col=cols)

## X axis: chromosome labels
axis(1, at=chr_mid, labels=fai\$CHR, tick=FALSE, cex.axis=0.8)
box()
__EOF__


mv TMP/${prefix}* ./


cat << EOF > versions.yml
${task.process}:
    R: \$(R --version | awk '(NR==1) {print \$3;}')
EOF
"""
stub:
    def format = task.ext.format?:"png"
    def prefix = task.ext.prefix?:"${assoc.baseName}" //not meta.id please
"""
touch versions.yml ${prefix}.${format}
"""
}
