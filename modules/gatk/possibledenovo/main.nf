/*

Copyright (c) 2025 Pierre Lindenbaum

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
process GATK_POSSIBLE_DENOVO {
    label "process_single"
    tag "${meta.id}"
    conda "${moduleDir}/../../../conda/bioinfo.02.yml"
    afterScript "rm -rf TMP"
    input:
        tuple val(meta1),path(fasta)
        tuple val(meta2),path(fai)
        tuple val(meta3),path(dict)
        tuple val(meta4),path(pedigree)
        tuple val(meta),path(vcf),path(vcfidx)
    output:
        tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
        path("versions.yml"),emit:versions
    script:
        def prefix = task.ext.prefix?:vcf.baseName+".denovo"
        def input_is_bcf = vcf.name.endsWith(".bcf")
    """
    mkdir -p TMP

    if ${input_is_bcf}
    then
        bcftools view -O z --threads ${task.cpus} -o TMP/jeter1.vcf.gz '${vcf}'
        bcftools index -f -t --threads ${task.cpus} TMP/jeter1.vcf.gz
    fi

    awk -f "${moduleDir}/pedigree4gatk.awk" "${pedigree}" > TMP/jeter.ped

    gatk --java-options "-XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" VariantAnnotator \\
        -R "${fasta}" \\
        --annotation PossibleDeNovo \\
        --pedigree  TMP/jeter.ped \\
        -V ${input_is_bcf? "TMP/jeter1.vcf.gz" : "\"${vcf}\""} \\
        -O "TMP/${prefix}.vcf.gz"

    rm -f TMP/jeter1.vcf.gz TMP/jeter1.vcf.gz.tbi

    bcftools index -f -t \\
        --threads "${task.cpus}" \\
        "TMP/${prefix}.vcf.gz"
    
    mv "TMP/${prefix}.vcf.gz" ./
    mv "TMP/${prefix}.vcf.gz.tbi" ./


cat << END_VERSIONS > versions.yml
"${task.process}":
	gatk: "\$(gatk --version 2>&1  | paste -s -d ' ' | tr -c -d 'A-Za-z0-9._-' )"
END_VERSIONS
    """
    }
