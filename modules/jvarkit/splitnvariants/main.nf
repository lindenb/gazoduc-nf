
process SPLIT_N_VARIANTS {
tag "${meta.id?:""}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta ),path(vcf),path(idx),path(optional_bed)
output:
	tuple val(meta),path("OUT/*.vcf.gz"),optional:true,emit:vcf
	tuple val(meta),path("OUT/*.vcf.gz.tbi"),optional:true,emit:tbi
	tuple val(meta),path("*.MF"),optional:true,emit:manifest
	path("versions.yml"),emit:version
script:
	def method = params.split_vcf_method?:""
	if(method.trim().isEmpty()) throw new IllegalArgumentException("method undefined for ${task.process}");
	def has_bed = optional_bed?true:false
	def args1 = task.ext.args1?:""
	def prefix = task.ext.prefix?:vcf.baseName+".split"
"""
hostname 1>&2
set -o pipefail
mkdir -p OUT TMP

if ${has_bed}
then
	bcftools index -s "${vcf}" |\\
		awk -F '\t' '{printf("%s\t0\t%s\\n",\$1,\$2);}' |\\
		LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter1.bed

		
	${optional_bed.endsWith(".gz")?"gunzip -c":"cat"} "${optional_bed}" |\\
		LC_ALL=C  sort -T TMP -t '\t' -k1,1 -k2,2n |\\
		bedtools merge > TMP/jeter2.bed

	bedtools intersect \\
		-a  TMP/jeter1.bed \\
		-b TMP/jeter2.bed > TMP/jeter3.bed

	# prevent empty bed
	if test ! -s TMP/jeter3.bed ; then
		echo -e 'K_\t0\t1' > TMP/jeter3.bed
	fi
fi

bcftools view ${args1} \\
	${has_bed?" --regions-file TMP/jeter3.bed":""} \\
	"${vcf}" |\\
	java -Xmx${task.memory.giga}g  -XX:-UsePerfData  -Djava.io.tmpdir=TMP vcfsplitnvariants \\
	--manifest ${prefix}.MF \\
	${method} \\
	-o OUT/${prefix}

find OUT/ -type f -name "*.vcf.gz" | while read F
do
	bcftools index --threads ${task.cpus} -t --force "\${F}"
done

cat << EOF > versions.yml
${task.process}:
	jvarkit: TODO
EOF
"""
}

