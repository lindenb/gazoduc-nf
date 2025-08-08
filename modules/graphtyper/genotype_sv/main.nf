
process GRAPHTYPER_GENOTYPE_SV {
label "process_single"
tag "${vcf.name}"
afterScript "rm -rf TMP sv_results"
conda "${moduleDir}/../../../conda/graphtyper.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
	tuple val(meta4),path("BAMS/*")
	tuple val(meta ),path(vcf),path(vcfidx)
output:
	tuple val(meta ), path("*.bcf"), path("*.csi"),emit:vcf
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:"\${MD5}"
"""
hostname 1>&2
mkdir -p TMP

export TMPDIR=\${PWD}/TMP

find BAMS/ \\( -name "*.cram" -o -name "*.bam" \\) | sort -T TMP -V > TMP/bams.list

bcftools query -f '%CHROM\t%POS0\t%END\\n' '${vcf}' |\\
    sort -T TMP -t '\t' -k1,1 -k2,2n |\\
    bedtools merge |\\
	awk -F '\t' '{printf("%s:%d-%s\\n",\$1,int(\$2)+1,\$3);}' > TMP/jeter.intervals

# prevent empty intervals
if test ! -s TMP/jeter.intervals
then
	awk '(NR==1){printf("%s:1-2\\n",\$1);}' "${fai}" > TMP/jeter.intervals
fi

MD5=`cat TMP/bams.list  TMP/jeter.intervals | md5sum | cut -d ' ' -f1`

graphtyper genotype_sv \\
	"${fasta}" \\
	"${vcf}" \\
	--force_no_copy_reference \\
	--force_use_input_ref_for_cram_reading \\
	--region_file TMP/jeter.intervals \\
	--sams=TMP/bams.list \\
	--threads=${task.cpus}

find \${PWD}/sv_results/ -type f -name "*.vcf.gz" | grep -v '/input_sites/' > TMP/vcf.list

bcftools concat \\
    --file-list TMP/vcf.list \\
	--allow-overlaps  \\
	--threads ${task.cpus} \\
    -O u |\\
	bcftools sort -T TMP/sort -O b -o TMP/jeter.bcf

bcftools index -f --threads ${task.cpus} TMP/jeter.bcf

mv  TMP/jeter.bcf  ./${prefix}.bcf
mv  TMP/jeter.bcf.csi  ./${prefix}.bcf.csi

cat << EOF > versions.yml
"${task.process}"
    graphtyper: todo
EOF
"""
}
