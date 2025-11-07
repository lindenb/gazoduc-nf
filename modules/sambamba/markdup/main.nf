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

process SAMBAMBA_MARKDUP {
tag "${meta.id?:""} ${bam.name}"
label "process_short"
conda "${moduleDir}/../../../conda/sambamba.yml"
afterScript "rm -rf TMP"

input:
	tuple val(meta),path(bam)
output:
	tuple val(meta),path("*.bam"),path("*.bai"),emit:bam
	path("versions.yml"),emit:versions
script:
	def overflow_size = task.ext.overflow_size?:"600000"
	def prefix = task.ext.prefix?:"${meta.id}.markdup"
	if(bam.name.endsWith(".cram")) throw new IllegalArgumentException("CRAM are not supported by sambamba");
"""
hostname 2>&1
mkdir TMP

# avoid too many files open. ulimit doesn't work...
ulimit -s unlimited || true

# markdup
sambamba markdup \\
	--overflow-list-size ${overflow_size} \\
	--tmpdir=TMP \\
	-t ${task.cpus} \\
	"${bam}" TMP/jeter2.bam
mv TMP/jeter2.bam TMP/jeter.bam

samtools index  -@ ${task.cpus} "TMP/jeter.bam"


mv TMP/jeter.bam "${prefix}.bam"
mv TMP/jeter.bam.bai "${prefix}.bam.bai"


cat << EOF > versions.yml
"${task.process}":
    sambamba: \$(sambamba --version 2>&1 | sort | uniq | paste -s -d ' ')
    samtools: \$(samtools  --version | head -n 1| cut -d ' ' -f2)
</properties>
EOF
"""
stub:
"""
touch  "${sample}.markdup.bam"   "${sample}.markdup.bam.bai"
echo "<properties/>" > version.xml
"""
}

