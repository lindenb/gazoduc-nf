process MERGE_GENOTYPES {
    tag "${meta.id}"
    label "process_short"
    afterScript "rm -rf TMP"
    conda "${moduleDir}/../../conda/delly.yml"
    input:
	    tuple val(meta),path("VCFS/*")
    output:
	    tuple val(meta),path("*.bcf",arity:"1"),path("*.csi",arity:"1"),emit:vcf
        path("versions.yml"),emit:versions
    script:
        def prefix = task.ext.prefix?:(meta.id?:"merged")
    // note to self , cd gazoduc-nf if not enough memory
    """
    hostname 1>&2
    ulimit -s unlimited || true
    mkdir -p TMP
    find VCFS/ -name "*.bcf" | sort > TMP/jeter.list
    test -s TMP/jeter.list

    bcftools merge --force-single --threads ${task.cpus} -m id -O b9 -o TMP/${prefix}.bcf --file-list TMP/jeter.list
    bcftools index --threads ${task.cpus} --csi TMP/${prefix}.bcf 
    
    mv TMP/${prefix}.bcf ./
    mv TMP/${prefix}.bcf.csi ./

cat << END_VERSIONS > versions.yml
"${task.process}":
    bcftools: \$(bcftools version | awk '(NR==1) {print \$NF;}')
END_VERSIONS
    """

stub:
 def prefix = task.ext.prefix?:(meta.id?:"merged")
"""
touch versions.yml {prefix}.bcf {prefix}.bcf.csi
"""
}