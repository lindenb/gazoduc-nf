
workflow  {
	if(params.vcf.endsWith(".list")) {
		ch2 = Channel.fromPath(params.vcf).splitText().map{file(it.trim())}
		}
	else
		{
		ch2 = Channel.fromPath(params.vcf)
		}
	ch3 = PLINK2_VCF2PGEN(ch2)
	ch4 = PLINK2_MERGE_PGEN(ch3.output.collect())
	ch5 = PLINK2_PGEN2BGEN(ch4.output)
	}

process PLINK2_VCF2PGEN {
label "process_short"
afterScript "rm -rf TMP"
tag "${vcf.name}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	path(vcf)
output:
	path("${vcf.baseName}.*"),emit:output
script:
	def args = " -set-id +'%VKX' "
"""
set -o pipefail
mkdir -p TMP

plink2 ${vcf.name.endsWith(".bcf")?"--bcf":"--vcf"} "${vcf}"  \\
	--double-id \\
	--make-pgen erase-phase \\
	--threads ${task.cpus} \\
	--out "${vcf.baseName}"
"""
}


process PLINK2_MERGE_PGEN {
label "process_short"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
        path("INPUT/*")
output:
        path("merged.*"),emit:output
script:
        def args = " -set-id +'%VKX' "
"""
set -o pipefail
mkdir -p TMP
find INPUT -type l -name "*.pgen" | sed 's/\\.pgen\$//' > TMP/jeter.list

if test `wc -l < TMP/jeter.list` -eq 1
then
	ln -s INPUT/*.log merged.log
	ln -s INPUT/*.pvar merged.pvar
	ln -s INPUT/*.psam merged.psam
	ln -s INPUT/*.pgen merged.pgen
else

plink2 --pmerge-list TMP/jeter.list \\
        --make-pgen \\
	--threads ${task.cpus} \\
        --out "merged"
fi
"""
}



process PLINK2_PGEN2BGEN {
label "process_short"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
        path("INPUT/*")
output:
        path("merged.*"),emit:output
script:
        def args = " -set-id +'%VKX' "
"""
set -o pipefail
mkdir -p TMP

plink2 \\
            --pfile INPUT/merged \\
            --export bgen-1.2 ref-first bits=8 id-paste=iid id-delim='-' \\
            --out "merged"

"""
}
