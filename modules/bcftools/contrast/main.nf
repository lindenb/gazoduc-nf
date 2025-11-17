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

process BCTOOLS_CONTRAST {
    label "process_single"
    tag "${meta.id}"
    conda "${moduleDir}/../../../conda/bioinfo.01.yml"
    afterScript "rm -rf TMP"
    input:
        tuple val(meta2),path(pedigree)
        tuple val(meta),val(vcf),path(vcfidx)
    output:
        tuple val(meta),path("*.bcf"),path("*.csi"),emit:vcf
        path("versions.yml"),emit:versions
    script:
        def prefix = task.ext.prefix?:vcf.baseName+".contrast"
    """
    mkdir -p TMP

bcftools query -l '${vcf}' | sort | uniq > TMP/jeter1.txt

awk '(\$6=="case" || \$6=="affected") {print \$2;}' '${pedigree}' | sort | uniq > TMP/jeter2.txt
comm -12 TMP/jeter1.txt  TMP/jeter2.txt > TMP/jeter.cases.txt

awk '(\$6=="control" || \$6=="unaffected") {print \$2;}' '${pedigree}' | sort | uniq > TMP/jeter2.txt
comm -12 TMP/jeter1.txt  TMP/jeter2.txt > TMP/jeter.ctrls.txt



if [  -s "TMP/jeter.cases.txt" ] && [  -s "TMP/jeter.ctrls.txt"	]
then

	bcftools +contrast \
		-0 TMP/jeter.ctrls.txt \
		-1 TMP/jeter.cases.txt \
		-a PASSOC,FASSOC,NASSOC,NOVELAL,NOVELGT -O b -o TMP/jeter.bcf '${vcf}'

    bcftools index --threads ${task.cpus} -f  TMP/jeter.bcf

else
     bcftools view --write-index --threads ${task.cpus} -o TMP/jeter.bcf '${vcf}'
fi

  
mv  TMP/jeter.bcf ${prefix}.bcf
mv  TMP/jeter.bcf.csi ${prefix}.bcf.csi


cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
stub:
     def prefix = task.ext.prefix?:vcf.baseName+".contrast"
"""
touch versions.yml ${prefix}.bcf ${prefix}.bcf.csi
"""
}
