include {k1_signature} from '../../../modules/utils/k1.nf'


process BCTOOLS_MENDELIAN2 {
    label "process_single"
    tag "${meta.id}"
    conda "${moduleDir}/../../../conda/bioinfo.01.yml"
    afterScript "rm -rf TMP"
    input:
        tuple val(meta1),path(fai)
        tuple val(meta2),path(pedigree)
        tuple val(meta),path(vcf),path(vcfidx)
    output:
        tuple val(meta),path("*{.vcf.gz,.bcf}"),path("*{.tbi,.csi}"),emit:vcf
        path("versions.yml"),emit:versions
    script:
        def k1 = k1_signature()
        def prefix = task.ext.prefix?:vcf.baseName+".mendelian"
        def suffix = task.ext.suffix?:".bcf"
	def test_trio_exists = task.ext.test_trio_exists?:true
    """
    
mkdir -p TMP


cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg38}\t--rules GRCh38
1:${k1.hg19}\t--rules GRCh37
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' |\\
    sed 's/^chr//' |\\
    sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv

join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv |\\
    sort |\\
    uniq > TMP/jeter.rule


    awk '{SEX=1;if(\$5=="female") SEX=2;printf("%s\t%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,\$4,SEX);}' ${pedigree} > TMP/jeter.ped
    cat TMP/jeter.ped 1>&2
    ${test_trio_exists ? "test -s TMP/jeter.ped" : ""}

    if test -s TMP/jeter.ped
    then


    bcftools +mendelian2 \\
        -O ${suffix.contains("bcf")?"b":"z"} \\
        -o "TMP/${prefix}${suffix}" \\
        `cat TMP/jeter.rule` \\
        -m a -P TMP/jeter.ped "${vcf}"

   else

	bcftools view  \\
		--threads "${task.cpus}" \\
		-O ${suffix.contains("bcf")?"b":"z"}  \\
		-o "TMP/${prefix}${suffix}"  "${vcf}"

   fi

    bcftools index \\
        -f \\
        ${suffix.contains("bcf")?"":"-t"} \\
        --threads "${task.cpus}" \\
        "TMP/${prefix}${suffix}"

    mv TMP/${prefix}* ./

cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
    """
    }
