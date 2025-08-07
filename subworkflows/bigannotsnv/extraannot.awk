/^#CHROM/ {
    printf("##INFO=<ID=MEDIAN_DP,Number=1,Type=Float,Description=\"Median DEPTH\">\n");
    printf("##INFO=<ID=MEAN_DP,Number=1,Type=Float,Description=\"Average DEPTH\">\n");
    printf("##INFO=<ID=MIN_DP,Number=1,Type=Float,Description=\"Min DEPTH\">\n");
    printf("##INFO=<ID=MAX_DP,Number=1,Type=Float,Description=\"Max DEPTH\">\n");


    printf("##INFO=<ID=MEDIAN_GQ,Number=1,Type=Float,Description=\"Median Genotype Quality\">\n");
    printf("##INFO=<ID=MEAN_GQ,Number=1,Type=Float,Description=\"Average Genotype Quality\">\n");
    printf("##INFO=<ID=MIN_GQ,Number=1,Type=Float,Description=\"Min Genotype Quality\">\n");
    printf("##INFO=<ID=MAX_GQ,Number=1,Type=Float,Description=\"Max Genotype Quality\">\n");

    printf("##INFO=<ID=MEDIAN_HET_AD,Number=1,Type=Float,Description=\"Median AD ratio for Heterozygous Genotypes (should be close to 0.5)\">\n");
    printf("##INFO=<ID=MEAN_HET_AD,Number=1,Type=Float,Description=\"Average AD ratio for Heterozygous Genotypes (should be close to 0.5)\">\n");
    

    printf("##INFO=<ID=MEDIAN_AA_AD,Number=1,Type=Float,Description=\"Median ratio between ALT and sum(AD) for HOM_VAR genotypes  (should be close to 1.0)\">\n");
    printf("##INFO=<ID=MEAN_AA_AD,Number=1,Type=Float,Description=\"Average ratio between ALT  and sum(AD) ralleles for HOM_VAR genotypes  (should be close to 1.0)\">\n");

    printf("##INFO=<ID=MEDIAN_RR_AD,Number=1,Type=Float,Description=\"Median ratio between REF and sum(AD) for HOM_REF genotypes  (should be close to 1.0)\">\n");
    printf("##INFO=<ID=MEAN_RR_AD,Number=1,Type=Float,Description=\"Average ratio between REF  and sum(AD) ralleles for HOM_REF genotypes  (should be close to 1.0)\">\n");

    printf("##INFO=<ID=PURE_HOM,Number=0,Type=Flag,Description=\"All Homozygous Genotypes have no trace of ALT allele\">\n");
    
    printf("##INFO=<ID=SINGLETON,Number=1,Type=String,Description=\"Name of the unique sample having ALT allele\">\n");
 
    printf("##INFO=<ID=N_ALT_SAMPLES,Number=1,Type=String,Description=\"Number of samples having a ALT allele\">\n");
    }

    {
    print;
    }