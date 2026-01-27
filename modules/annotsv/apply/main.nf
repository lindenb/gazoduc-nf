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
/**


COMMAND LINE USAGE

       $ANNOTSV/bin/AnnotSV -SVinputFile 'Path of your VCF or BED input file with SV coordinates' >& AnnotSV.log &


OPTIONS

-annotationsDir:                Path of the annotations directory

-annotationMode:                Description of the types of lines produced by AnnotSV
                                Values: both (default), full or split

-bcftools:                      Path of the bcftools local installation

-bedtools:                      Path of the bedtools local installation

-benignAF:                      Allele frequency threshold to select the benign SV in the data sources
                                Range values: [0.001-0.1], default = 0.01 (i.e. 1%)	

-candidateGenesFile:            Path of a file containing the candidate genes of the user (gene names can be space-separated, tabulation-separated, or line-break-separated)

-candidateGenesFiltering:       To select only the SV annotations ("split" and "full") overlapping a gene from the "candidateGenesFile"
                                Values: 0 (default) or 1

-candidateSnvIndelFiles:        Path of the filtered VCF input file(s) with SNV/indel coordinates for compound heterozygotes report (optional)
                                Gzipped VCF files are supported as well as regular expression

-candidateSnvIndelSamples:      To specifiy the sample names from the VCF files defined with the -candidateSnvIndelFiles option (sample names can be coma-separated or semicolon -separated)
                                Default: use all samples from the filtered VCF files

-genomeBuild:                   Genome build used
                                Values: GRCh38 (default) or GRCh37 or CHM13 or or mm9 or mm10

-help:                          More information on the arguments

-hpo:                           HPO terms list describing the phenotype of the individual being investigated
                                Values: use comma, semicolon or space separated class values
                                Default = "" (e.g.: "HP:0001156,HP:0001363,HP:0011304")

-includeCI:                     To expand the "start" and "end" SV positions with the VCF confidence intervals (CIPOS, CIEND) around the breakpoints
                                AnnotSV keeps the CIPOS and CIEND information that comes first in the INFO column (even if the fields are CIPOS95, CIEND95 or tool_CIPOS, tool_CIEND).
                                Values: 1 (default) or 0

-metrics:                       Changing numerical values from frequencies to us or fr metrics (e.g. 0.2 or 0,2)
                                Values: us (default) or fr

-minTotalNumber:                Minimum number of individuals tested to consider a benign SV for the ranking
                                Range values: [100-1000], default = 500

-missingGTinSamplesid:          To report sample IDs with missing alleles in the GT field (./. or .|.) in the "Samples_ID" output field
                                Values: 1 (default) or 0

-outputDir:                     Output path name

-outputFile:                    Output path and file name 

-overlap:                       Minimum overlap (%) between user features (User BED) and the annotated SV to be reported
                                Range values: [0-100], default = 100

-overwrite:                     To overwrite existing output results
                                Values: 1 (default) or 0

-promoterSize:                  Number of bases upstream from the transcription start site
                                Default = 500

-rankFiltering:                 To select the SV of a user-defined specific class (from 1 to 5; or NA)
                                Values: use comma separated class values, or use a dash to denote a range of values
                                (e.g.: "3,4,5" or "3-5"), default = "1-5,NA"

-reciprocal:                    Use of a reciprocal overlap between SV and user features (only for annotations with features overlapping the SV)
                                Values: 0 (default) or 1

-REreport:                      Create a report to link the annotated SV and the overlapped regulatory elements (coordinates and sources)
                                Values: 0 (default) or 1

-REselect1:                     To report only the morbid, HI, TS, candidate and phenotype matched genes 
                                Values: 1 (default) or 0
	
-REselect2:                     To report only the genes not present in "Gene_name"
                                Values: 1 (default) or 0

-samplesidBEDcol:               Number of the column reporting the samples ID for which the SV was called (if the input SV file is a BED)
                               	Range values: [4-[, default = -1 (value not given)
                                (Samples ID should be comma or space separated) 

-snvIndelFiles:	                Path of the VCF input file(s) with SNV/indel coordinates used for false positive discovery
                                Use counts of the homozygous and heterozygous variants
                                Gzipped VCF files are supported as well as regular expression

-snvIndelPASS:                  Boolean. To only use variants from VCF input files that passed all filters during the calling (FILTER column value equal to PASS)
                                Values: 0 (default) or 1

-snvIndelSamples:               To specify the sample names from the VCF files defined from the -snvIndelFiles option
                                Default: use all samples from the VCF files

-SVinputFile:                   Path of the input file (VCF or BED) with SV coordinates
                                Gzipped VCF file is supported

-SVinputInfo:                   To extract the additional SV input fields and insert the data in the outputfile
                                Values: 1 (default) or 0

-SVminSize:                     SV minimum size (in bp)
                                AnnotSV does not annotate small deletion, insertion and duplication from a VCF input file 
                                Default = 50

-svtBEDcol:                     Number of the column describing the SV type (DEL, DUP) if the input SV file is a BED
                                Range values: [4-[, default = -1 (value not given)

-tx:                            Origin of the transcripts (RefSeq or ENSEMBL)
                                Values: RefSeq (default) or ENSEMBL

-txFile:                        Path of a file containing a list of preferred genes transcripts to be used in priority during the annotation (Preferred genes transcripts names should be tab or space separated)

-variantconvertDir:             Path of the variantconvert directory
                                (by default, the variantconvert tool distributed by annotSV is used)

-variantconvertMode:            variantconvert conversion mode
                                Values: combined (default) or full or fullsplit

-version:                       Version of the AnnotSV program

-vcf:                           Creation of a VCF output file format (-svtBEDcol needs to be defined too)
                                Values: 0 (default) or 1



*/
process ANNOTSV_APPLY {
tag "${meta.id?:""}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/annotsv.yml"
input:
    tuple val(meta1),path(annotsv_db)
    tuple val(meta2),path(dict) //for build/ucsc_name
    tuple val(meta),path(vcf)
output:
	tuple val(meta),path("*.annotSV.tsv"),optional:true,emit:annot
        tuple val(meta),path("*.unannotated.tsv"),optional:true,emit:unannotated
	path("versions.yml"),emit:versions
script:
   def prefix = task.ext.prefix?:"${meta.id}.annotsv"
   def ucsc_name = task.ext.ucsc_name?:"${meta2.ucsc_name}"
   def build = task.ext.build?:(ucsc_name=="hg38"?"GRCh38":(ucsc_name=="hg19"?"GRCh37":"unknown"))
   def args1 = task.ext.args1?:""
"""
mkdir -p TMP

# run annotSV
AnnotSV \\
	-annotationsDir "\${PWD}/${annotsv_db.name}" \\
	${args1} \\
	-genomeBuild ${build} \\
	-SVinputFile "${vcf}" \\
	-outputDir "TMP" \\
	-outputFile "${prefix}.annotSV.tsv" 1>&2

find TMP -type f 1>&2
mv -v TMP/*.tsv ./

touch versions.yml
"""

stub:
     def prefix = task.ext.prefix?:"${meta.id}.annotsv" 
"""
mkdir -p OUT
touch versions.yml OUT/${prefix}.annotSV.tsv OUT/${prefix}.annotSV.unannotated.tsv 
"""
}
