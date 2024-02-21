/*

Copyright (c) 2024 Pierre Lindenbaum

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
		burden_ch = BURDEN_SETFILE(params.genomeId, file(params.vcf), file(params.pedigree),file(params.bed))
		}

workflow BURDEN_SETFILE {
	take:
		genomeId
		vcf
		pedigree
		bed
	main:

		version_ch = Channel.empty()
		to_zip = Channel.empty()



		header_ch = RVTESTS_REHEADER_01([:], genomeId)
		version_ch = version_ch.mix(header_ch.version)
		ped2_ch = PEDIGREE_FOR_RVTESTS([:],pedigree)


		vcf2bed_ch = VCF_TO_BED([:],vcf)

		bed_ch = DIGEST_BED(genomeId,bed,vcf2bed_ch.bed)
	



		assoc_ch = RVTEST_BY_ID(
			genomeId,
			vcf2bed_ch.bed,
			header_ch.output,
			ped2_ch.pedigree,
			bed_ch.bed,
			bed_ch.interval_ids.splitText().map{it.trim()}
			)
		//version_ch = version_ch.mix(assoc_ch.version)
		
		concat_ch = CONCAT_FILES_01([:],assoc_ch.output.collect())
		version_ch = version_ch.mix(concat_ch.version)

		digest_ch = RVTESTS_POST_PROCESS([:], genomeId, "${params.prefix?:""}combinations" ,concat_ch.output)
                version_ch = version_ch.mix(digest_ch.version)


                params4multiqc_ch = PARAMS_MULTIQC([:])

                multiqc_ch = MULTIQC_01(["title":"${params.prefix}Combinations","comment":"VCF ${params.vcf}"], digest_ch.to_multiqc.concat(params4multiqc_ch.output).collect());
                version_ch = version_ch.mix(multiqc_ch.version)

                version_ch = MERGE_VERSION("burden coding",version_ch.collect())
	}


process DIGEST_BED {
tag "${bed.name}"
afterScript "rm -rf TMP"
memory "2G"
input:
	val(genomeId)
	path(bed)
	path(vcf2bed)
output:
	path("chromosomes.txt"),emit:chromosomes
	path("ids.txt"),emit:interval_ids
	path("interval_set.bed"),emit:bed
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
"""
hostname 1>&2
${moduleLoad("bedtools jvarkit")}
mkdir -p TMP
set -x

${bed.name.endsWith(".gz")?"gunzip -c":"cat"} "${bed}" |\
	tr -d '\\r' |\
	grep -v "^#" |\
	cut -f1,2,3,4 |\
	java  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -R "${reference}" --column 1 --convert SKIP |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\
	uniq > TMP/jeter.bed

# check output is not empty
test -s TMP/jeter.bed


# VCF bed
cut -f1,2,3 "${vcf2bed}" |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\
	uniq > TMP/jeter2.bed


bedtools intersect -a TMP/jeter.bed -b TMP/jeter2.bed -u  > TMP/jeter3.bed
mv -v TMP/jeter3.bed TMP/jeter.bed
echo "the following test fails if there is no overlap between the regions of interest and the VCF(s)" 1>&2
test -s TMP/jeter.bed

# merged, currently not used...
cut -f1,2,3 TMP/jeter.bed | bedtools merge > TMP/merged.bed

# check fourth column is ok
awk -F '\t' '!(\$4 ~ /^[A-Z\\-\\.0-9a-z_]+\$/)' TMP/jeter.bed |head | cat > TMP/jeter2.bed
echo "The following test will fail if the fourth column has a bad syntax" 1>&2
head TMP/jeter2.bed
test ! -s TMP/jeter2.bed

cut -f 1 TMP/jeter.bed | uniq | sort | uniq > TMP/chromosomes.txt

cut -f 4 TMP/jeter.bed | uniq | sort |	uniq > TMP/ids.txt



mv -v TMP/chromosomes.txt ./
mv -v TMP/ids.txt ./
mv -v TMP/jeter.bed interval_set.bed
"""
}

process RVTEST_BY_ID {
tag "${intervalid}"
afterScript "rm -rf TMP"
input:
	val(genomeId)
	path(vcf2bed)
	path(reheader)
	path(pedigree)
	path(bed)
	val(intervalid)
output:
	path("assoc.list"),emit:output
when:
	true
script:
	def  rvtest_params = params.rvtests.arguments
	def winsize = params.winsize
	def winshift = params.winshift
	def with_makewindows = winsize>0 && winshift>0
"""
hostname 1>&2
${moduleLoad("rvtests bcftools bedtools")}
set -o pipefail
set -x

mkdir -p ASSOC TMP


awk -F '\t' '(\$4=="${intervalid}")' "${bed}" |\
	cut -f1,2,3 |\
	sort -T TMP -t '\t' -k1,1 -k2,2n |\
	bedtools merge > TMP/intervals.01.bed

test -s TMP/intervals.01.bed

# only split record if the cumulative length is greater than winsize
if ${with_makewindows} && awk -F '\t' 'BEGIN{N=0;} {N+=int(\$3)-int(\$2);} END {if(N>${winsize}) print "LARGE";}' TMP/intervals.01.bed | grep LARGE 1>&2  ; then


cat << __EOF__ > TMP/jeter.py
import sys

class StartEnd:
    def __init__(self,contig,start,end):
        self.contig = contig
        self.start = start
        self.end = end

def read_bed_file_from_stdin():
    start_end_list = []
    for line in sys.stdin:
        if line.startswith('#'):
            continue  # Ignore comment lines
        parts = line.strip().split('\\t')
        if len(parts) < 3:
            continue  # Skip malformed lines
        try:
            start = int(parts[1])
            end = int(parts[2])
            start_end_list.append(StartEnd(parts[0],start, end))
        except ValueError:
            continue  # Skip lines where start/end are not integers
    start_end_list.sort(key=lambda x: x.start)
    return start_end_list

label="${intervalid}"
setfile=[]
winsize=${winsize}
winshift=${winshift}
L = read_bed_file_from_stdin()
start = L[0].start
while start < L[len(L)-1].end:
    last = start
    remain=winsize
    i=0
    s=[]
    while remain>0 and i < len(L):
        if L[i].end < start:
            i=i+1
            continue
        x1 = max(L[i].start,start)
        x2 = min(L[i].end,x1+remain)
        contig = L[0].contig
        if contig.startswith("chr"):
            contig = contig[3:]
        s.append(contig+":"+str(x1+1)+"-"+str(x2))
        i=i+1
        last = x2
        remain -= x2-x1
    setfile.append(str(start)+"_"+str(last)+"\t"+ (",".join(s)))
    start+=winshift

if len(setfile)==1:
    print(label+"_"+setfile[0])
else:
    for i,s in enumerate(setfile):
        print(label+"_"+s)
__EOF__


	python3 TMP/jeter.py < TMP/intervals.01.bed > TMP/jeter.setfile

else

	#create setfile with normalized chr
	echo -n '${intervalid}\t' > TMP/jeter.setfile
	sed 's/^chr//' TMP/intervals.01.bed | awk -F '\t' '{printf("%s:%d-%d\\n",\$1,int(\$2)+1,\$3);}' | paste -sd, >> TMP/jeter.setfile

fi

#find vcf(s) for that intervals
sort -T TMP -t '\t' -k1,1 -k2,2n "${vcf2bed}" |\\
	bedtools intersect -a - -b TMP/intervals.01.bed -u |\\
	cut -f4 | sort | uniq > TMP/vcfs.list

# no vcf for those data. use a random
if test ! -s TMP/vcfs.list ; then
	head -n1 "${vcf2bed}" | cut -f 4 > TMP/vcfs.list
fi

# paranoid
test -s TMP/vcfs.list

# generate VCF for the setfile
bcftools concat -a --regions-file TMP/intervals.01.bed --file-list TMP/vcfs.list -O u |\
        bcftools annotate  --rename-chrs "${reheader}" -O z -o TMP/jeter.vcf.gz

bcftools index --force -t TMP/jeter.vcf.gz	

rvtest  --noweb \
        --inVcf TMP/jeter.vcf.gz \
	--setFile "TMP/jeter.setfile" \
	--pheno "${pedigree}" \
	--out "ASSOC/part" \
	${rvtest_params} 1>&2 2> TMP/last.rvtest.log


find \${PWD}/ASSOC -type f -name "part*.assoc" >  assoc.list
"""
}




runOnComplete(workflow);
