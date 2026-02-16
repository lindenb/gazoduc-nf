# PIHAT workflow

The probable relatives and duplicate samples are detected based on pairwise identity-by-state (IBS) from which a variable called PIHAT is calculated via PLINK. A higher PIHAT pairwise score indicates that the individuals are genetically more closely related.  
Several plots are also outputted: PCA and boxplot of average PIHAT scores per sample.

# Parameters

| **Option** | Status | Value |
|---|---|---|
|`--vcf (file)` | Required | A file with the `.list` suffix containing the path to your `vcf.gz` (multiple paths not accepted). Or it can be a CSV file with the following header: `vcf,tbi`. Or just the path to file with the `.vcf.gz` suffix.|
|`--fasta (file)` | Required | Fasta reference file associated to your vcf file.|
|`--onekg (file)` | Optional | A file with the `.list` suffix containing the path to the 1000 Genomes Project `vcf.gz` (multiple paths accepted). Or it can be a CSV file with the following header: `vcf,tbi`. Or just the path to file with the `.vcf.gz` suffix.|
|`--gnomad (file)` | Optional | A file with the `.list` suffix containing the path to the gnomad `vcf.gz`. Or it can be a CSV file with the following header: `vcf,tbi`. Or just the path to file with the `.vcf.gz` suffix.|
|`--sample2pop (file)` | Optional | TSV file mapping sample to population name (format: `sample_name  population`).|
|`--remove_samples (file)` | Optional | Text file containing samples names to remove from your vcf, one per line.|
|`--exclude_bed (file)` | Optional | BED file of loci to exclude from the analysis.|
