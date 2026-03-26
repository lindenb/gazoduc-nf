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

/*
>>> 2
$1        DATE ADDED TO CATALOG : 2008-06-16
$2                     PUBMEDID : 17434096
$3                 FIRST AUTHOR : Matarin M
$4                         DATE : 2007-05-06
$5                      JOURNAL : Lancet Neurol
$6                         LINK : www.ncbi.nlm.nih.gov/pubmed/17434096
$7                        STUDY : A genome-wide genotyping study in patients with ischaemic stroke: initial analysis and data release.
$8                DISEASE/TRAIT : Stroke
$9          INITIAL SAMPLE SIZE : 249 European ancestry cases, 268 European ancestry controls
$10     REPLICATION SAMPLE SIZE : NA
$11                      REGION : 18p11.21
$12                      CHR_ID : 18
$13                     CHR_POS : 11987273
$14            REPORTED GENE(S) : IMPA2
$15                 MAPPED_GENE : IMPA2
$16            UPSTREAM_GENE_ID : 
$17          DOWNSTREAM_GENE_ID : 
$18                SNP_GENE_IDS : ENSG00000141401
$19      UPSTREAM_GENE_DISTANCE : 
$20    DOWNSTREAM_GENE_DISTANCE : 
$21   STRONGEST SNP-RISK ALLELE : rs7506045-?
$22                        SNPS : rs7506045
$23                      MERGED : 0
$24              SNP_ID_CURRENT : 7506045
$25                     CONTEXT : intron_variant
$26                  INTERGENIC : 0
$27       RISK ALLELE FREQUENCY : 0.10
$28                     P-VALUE : 7E-7
$29                 PVALUE_MLOG : 6.154901959985743
$30              P-VALUE (TEXT) : 
$31                  OR or BETA : 5.39
$32               95% CI (TEXT) : [2.77-10.5]
$33  PLATFORM [SNPS PASSING QC] : Illumina [408803]
$34                         CNV : N
$35                MAPPED_TRAIT : stroke
$36            MAPPED_TRAIT_URI : http://www.ebi.ac.uk/efo/EFO_0000712
$37             STUDY ACCESSION : GCST000032
 : Genome-wide genotyping array
<<< 2

*/

process GWAS_CATALOG_DOWNLOAD_ASSOC{
tag "${meta.id}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
val(meta)
output:
	tuple val(meta),path("*.tsv"),emit:catalog
	path("versions.yml"),emit:versions
script:
	def url=task.ext.url?:"https://www.ebi.ac.uk/gwas/api/search/downloads/associations/v1.0.2?split=false"
	def prefix = task.ext.prefix?:"gwas-catalog-download-associations-alt-full"
"""
hostname 1>&2
mkdir -p TMP


curl -L -o TMP/jeter.zip "${url}"
unzip -p TMP/jeter.zip gwas-catalog-download-associations-alt-full.tsv > TMP/gwas-catalog-download-associations-alt-full.tsv


mv TMP/gwas-catalog-download-associations-alt-full.tsv ./${prefix}.tsv

cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${url}"
END_VERSIONS
"""

stub:
	def prefix = task.ext.prefix?:"${meta.id}.gwascatalog"
"""
touch versions.yml ${prefix}.tsv
"""
}
