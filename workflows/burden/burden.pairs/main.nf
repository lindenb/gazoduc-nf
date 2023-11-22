/*

Copyright (c) 2023 Pierre Lindenbaum

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
nextflow.enable.dsl=2
include {dumpParams;runOnComplete;moduleLoad;getVersionCmd} from '../../../modules/utils/functions.nf'
include {VCF_TO_BED} from '../../../modules/bcftools/vcf2bed.01.nf'
include {PEDIGREE_FOR_RVTESTS} from '../../../modules/rvtests/rvtests.cases.controls.ped.01.nf'
include {RVTESTS_REHEADER_01} from '../../../modules/rvtests/rvtests.reheader.01.nf'
include {RVTESTS_POST_PROCESS} from '../../../subworkflows/rvtests/rvtests.post.process.01.nf'
include {CONCAT_FILES_01} from '../../../modules/utils/concat.files.nf'
include {MULTIQC_01} from '../../../modules/multiqc/multiqc.01.nf'
include {PARAMS_MULTIQC} from '../../../modules/utils/params.multiqc.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.02.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'


if(params.help) {
        dumpParams(params);
        exit 0
        }
else
        {
        dumpParams(params);
        }



workflow {
		burden_ch = BURDEN_PAIRS([:], params.genomeId, file(params.vcf), file(params.pedigree),Channel.fromPath(params.setfile))
		}

workflow BURDEN_PAIRS {
	take:
		meta
		genomeId
		vcf
		pedigree
		setfile
	main:

		version_ch = Channel.empty()
		to_zip = Channel.empty()
		setrecords = setfile.splitCsv(sep:'\t',header:false).map{T->[name:T[0],intervals:T[1]]}

		def each_level = String.valueOf(params.levels).split("[,]+");
		combinations_ch = Channel.empty()		

		if(each_level.contains("1")) {
			combinations_ch = combinations_ch.mix(setrecords.map{T->[[T.name],T.intervals]})
			}
		if(each_level.contains("2")) {
			combinations_ch = combinations_ch.mix(
				setrecords.combine( setrecords ).
					filter{T->T[0].name.compareTo(T[1].name)<0}.
					map{T->[ [T[0].name,T[1].name],T[0].intervals+","+T[1].intervals ]}
				)
				
			}
		if(each_level.contains("3")) {
			combinations_ch = combinations_ch.mix(
				setrecords.combine( setrecords ).
					filter{T->T[0].name.compareTo(T[1].name)<0}.
					combine( setrecords ).
					filter{T->T[1].name.compareTo(T[2].name)<0}.
					map{T->[ [T[0].name,T[1].name,T[2].name],T[0].intervals+","+T[1].intervals+","+T[2].intervals ]}
				)
			}
		


		header_ch = RVTESTS_REHEADER_01([:], genomeId)
		version_ch = version_ch.mix(header_ch.version)
		ped2_ch = PEDIGREE_FOR_RVTESTS([:],pedigree)


		vcf2bed_ch = VCF_TO_BED([:],vcf)


		assoc_ch = RVTEST_BY_COMBINATION(
			[:],
			genomeId,
			vcf2bed_ch.bed,
			header_ch.output,
			ped2_ch.pedigree,
			combinations_ch
			)
		version_ch = version_ch.mix(assoc_ch.version)
		
		concat_ch = CONCAT_FILES_01([:],assoc_ch.output.collect())
		version_ch = version_ch.mix(concat_ch.version)

		digest_ch = RVTESTS_POST_PROCESS([:], genomeId, "${params.prefix?:""}combinations" ,concat_ch.output)
                version_ch = version_ch.mix(digest_ch.version)


                params4multiqc_ch = PARAMS_MULTIQC([:])

                multiqc_ch = MULTIQC_01(["title":"${params.prefix}Combinations","comment":"VCF ${params.vcf}"], digest_ch.to_multiqc.concat(params4multiqc_ch.output).collect());
                version_ch = version_ch.mix(multiqc_ch.version)

                version_ch = MERGE_VERSION("burden coding",version_ch.collect())

                html = VERSION_TO_HTML(version_ch.version)
                to_zip = to_zip.mix(html.html)

		
		/*
		to_zip = to_zip.mix(digest_ch.zip)

		version_ch = MERGE_VERSION(meta, "burden 1st intron", "Burden 1st intron ${vcf}", version_ch.collect())
		to_zip = to_zip.mix(version_ch.version)

		html = VERSION_TO_HTML(params,version_ch.version)
		to_zip = to_zip.mix(html.html)
		
	emit:
		version = version_ch
		zip = to_zip*/
	}


process RVTEST_BY_COMBINATION {
tag "${genes.join(" / ")} N=${genes.size()}"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(genomeId)
	path(vcf2bed)
	path(reheader)
	path(pedigree)
	tuple val(genes),val(intervals)
output:
	path("assoc.list"),emit:output
	path("version.xml"),emit:version
when:
	true
script:
	def  rvtest_params = params.rvtests.arguments
"""
hostname 1>&2
${moduleLoad("rvtests bcftools bedtools")}
set -o pipefail
set -x
mkdir -p ASSOC TMP


# convert comma separated setfile to BED
echo '${intervals}' |\
	tr "," "\\n" | tr ":-" "\t" |\
	awk -F '\t' '(NF==3) {printf("%s\t%d\t%d\\n",\$1,int(\$2)-1,\$3);}' |\
	sort -T TMP -t '\t' -k1,1 -k2,2n |\
	bedtools merge > TMP/intervals.bed
test -s TMP/intervals.bed


#create setfile with normalized chr
echo -n '${genes.join("_")}\t' > TMP/jeter.setfile
sed 's/^chr//' TMP/intervals.bed | awk -F '\t' '{printf("%s:%d-%d\\n",\$1,int(\$2)+1,\$3);}' | paste -sd, >> TMP/jeter.setfile

#find vcf(s) for that intervals
sort -T TMP -t '\t' -k1,1 -k2,2n "${vcf2bed}" |\\
	bedtools intersect -a - -b TMP/intervals.bed -u |\\
	cut -f4 | sort | uniq > TMP/vcfs.list

# no vcf for those data. use a random
if test ! -s TMP/vcfs.list ; then
	head -n1 "${vcf2bed}" | cut -f 4 > TMP/vcfs.list
fi

# paranoid
test -s TMP/vcfs.list

# generate VCF for the setfile
bcftools concat -a --regions-file TMP/intervals.bed --file-list TMP/vcfs.list -O u |\
        bcftools annotate  --rename-chrs "${reheader}" -O z -o TMP/jeter.vcf.gz
bcftools index -t TMP/jeter.vcf.gz
	

rvtest  --noweb \
        --inVcf TMP/jeter.vcf.gz \
	--setFile "TMP/jeter.setfile" \
	--pheno "${pedigree}" \
	--out "ASSOC/part" \
	${rvtest_params} 1>&2 2> TMP/last.rvtest.log


find \${PWD}/ASSOC -type f -name "part*.assoc" >  assoc.list


###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">rvtest for set of genes</entry>
	<entry key="rvtest.path">\$(which rvtest)</entry>
	<entry key="versions">${getVersionCmd("rvtest")}</entry>
</properties>
EOF
"""
}




runOnComplete(workflow);
