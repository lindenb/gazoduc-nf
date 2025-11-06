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
//  merge many files: https://github.com/dellytools/delly/issues/158
process MERGE_DELLY {
    label "process_short"
    afterScript "rm -rf TMP"
    conda "${moduleDir}/../../conda/delly2.yml"
    input:
        tuple val(meta1),path(fasta)
        tuple val(meta2),path(fai)
        tuple val(meta ),path("VCFS/*")
    output:
	     tuple val(meta ),path("*.bcf"),path("*.csi"),emit:vcf
         path("versions.yml"),emit:versions
    script:
        def prefix = task.prefix?:meta.id
	def with_bnd = (task.ext.with_bnd?:true).toBoolean()
    """
    hostname 1>&2
    export LC_ALL=C
    mkdir -p TMP
    export TMPDIR=\${PWD}/TMP
    export PATH=\${PWD}:\${PATH}

    # see https://github.com/dellytools/delly/issues/158

    find VCFS -name "*.bcf" > TMP/jeter.tsv
    test -s TMP/jeter.tsv    

    delly merge -o TMP/merged.bcf TMP/jeter.tsv 1>&2

	if ${!with_bnd} ; then
        bcftools view --threads ${task.cpus} -e 'INFO/SVTYPE="BND"' -O b -o ./${prefix}.bcf TMP/merged.bcf
		bcftools index --threads ${task.cpus}  -f ./${prefix}.bcf
	else
		mv TMP/merged.bcf ./${prefix}.bcf
		mv TMP/merged.bcf.csi ./${prefix}.bcf.csi
	fi

    touch versions.yml
    """

    stub:
    def prefix = task.prefix?:meta.id
    """
    touch versions.yml ${prefix}.bcf ${prefix}.bcf.csi
    """
    }
