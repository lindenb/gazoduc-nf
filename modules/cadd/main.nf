process CADD {
tag "${meta.id?:vcf.name}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"

input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta4),path(cadd_tabix),path(cadd_tbi)
    tuple val(meta),path(vcf),path(idx)
output:
    tuple val(meta),path("*.bcf"),path("*.csi"),emit:vcf
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:vcf.baseName+".cadd"
    def TAG = task.ext.tag?:"CADD"
"""
hostname 1>&2
mkdir -p TMP OUTPUT
set -o pipefail

export TMPDIR=\${PWD}/TMP

# get intervals for this VCF, remove chr prefix, sort and merge
bcftools query \\
    -f '%CHROM\t%POS0\t%END\\n' '${vcf}' |\\
    sed 's/^chr//' |\\
    LC_ALL=C sort -T TMP -t '\t' -S ${task.memory.kilo} -k1,1 -k2,2n |\\
    bedtools merge > TMP/intervals.bed

# check not empty
if test ! -s TMP/intervals.bed
then
    echo "chr22\t0\t1" > TMP/intervals.bed
fi

# create a SED file to convert the chromosomes
cut -f1 TMP/intervals.bed |\\
    sort | uniq |\\
    awk -F '\t' '{printf("%s\t%s\\n",\$1,\$1);}' |\\
    jvarkit  -Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP bedrenamechr -f "${fasta}" --column 2 --convert SKIP  |\\
    awk -F '\t' '{printf("s|^%s\t|%s\t|\\n",\$1,\$2);}' > TMP/jeter.sed
    
# check not empty
if test ! -s  TMP/jeter.sed
then
    echo "s/^ / /" > TMP/intervals.bed
fi

tabix -h  \\
    --regions TMP/intervals.bed \\
    "${cadd_tabix}" |\\
    sed -f TMP/jeter.sed |\\
    bgzip > TMP/database.tsv.gz

tabix -c 1 -c '#' -b 2 -e 2 TMP/database.tsv.gz

rm TMP/intervals.bed

bcftools annotate -c 'CHROM,POS,REF,ALT,${TAG}_RAWSCORE,${TAG}_PHRED' \\
    --threads ${task.cpus} \\
    -a TMP/database.tsv.gz \\
    -H '##INFO=<ID=${TAG}_RAWSCORE,Number=1,Type=Float,Description="Raw Score in CADD ${cadd_tabix}">' \\
    -H '##INFO=<ID=${TAG}_PHRED,Number=1,Type=Float,Description="Phred Score in CADD ${cadd_tabix}">' \\
    -O b -o  TMP/${prefix}.bcf \\
    '${vcf}'

bcftools index \\
    --threads ${task.cpus} \\
    --force TMP/${prefix}.bcf

mv TMP/${prefix}.bcf ./
mv TMP/${prefix}.bcf.csi ./

cat << END_VERSIONS > versions.yml
"${task.process}":
	TODO: "TODO"
END_VERSIONS
"""
}
