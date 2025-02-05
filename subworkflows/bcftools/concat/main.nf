def ALL_CONTIGS="ALL"

List makeSQRT(def L1) {
        def key = L1.get(0);
        def L = L1.get(1);
        int n = (int)Math.ceil(Math.sqrt(L.size()));
        if(n<25) n=25;
        def returnList = [];
        def currList = [];
        int i=0;
        for(;;) {
                if(i<L.size()) currList.add(L.get(i));
                if(i==L.size() || currList.size()==n) {
                        if(!currList.isEmpty()) returnList.add([key,currList]);
                        if(i==L.size()) break;
                        currList=[];
                        }
                i++;
                }
        return returnList;
        }

workflow BCFTOOLS_CONCAT {
take:
	vcfs // tuple [vcf,idx]
	bed //optional bed file
main:
	ch1 = CONTIGS(vcfs,bed)
	ch2 = ch1.output.splitText().
		map{[it[0].trim(),[it[1],it[2]]]}.
		flatMap{makeSQRT(it)}.
		map{[it[0],it[1].flatten()]}
	ch3 = LEVEL1(ch2,bed)
	ch4 = LEVEL2(ch3.output.groupTuple().map{[it[0],it[1].flatten()]})

emit:	
	output = ch4.output

}

process CONTIGS {
tag "${vcf.name}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
label "process_single"
input:
	tuple path(vcf),path(idx)
	path(bed)
output:
	tuple path("chroms.txt"),path(vcf),path(idx),emit:output
script:
	def method = task.ext.method?:"per_contig"

if(method.equals("per_contig") && !bed.name.contains("."))
"""
set -o pipefail
bcftools index -s "${vcf}" | cut -f1 > "chroms.txt"
"""
else if(method.equals("per_contig"))
"""
set -o pipefail
cut -f1 '${bed}' | uniq | sort | uniq > a.txt
bcftools index -s "${vcf}" | cut -f1 | sort | uniq > b.txt
comm -12 a.txt b.txt > "chroms.txt"
rm a.txt b.txt
"""
else if(method.equals("all"))
"""
echo "${ALL_CONTIGS}" > chroms.txt
"""
else
"""
echo "not a valid method ${method}"
"""
}


process LEVEL1 {
tag "${ctg} ${bed.name}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
label "process_short"
afterScript "rm -rf TMP"
input:
        tuple val(ctg),path("VCF/*")
	path(bed)
output:
	tuple val(ctg),path("concat.*"),emit:output
script:
	def args = "-a --rm-dups all"
"""
        mkdir -p TMP
	find VCF -type l \\( -name "*.bcf" -o -name "*.vcf.gz" \\)  | sort > TMP/jeter.list

	if ${!ctg.equals(ALL_CONTIGS) && bed.name.contains(".")}
	then
		awk -F '\t' '(\$1=="${ctg}") {print;N++;} END{if(N==0) printf("chrxxxxxx\t0\t1\\n");}' '${bed}' > TMP/jeter.bed

	elif ${bed.name.contains(".")}
	then
		cp '${bed}' TMP/jeter.bed	
	fi


        bcftools concat --threads ${task.cpus} \\
		${args} \\
		${!bed.name.contains(".")?(ctg.equals(ALL_CONTIGS)?"":"--regions \"${ctg}\""):"--regions-file TMP/jeter.bed"} \\
                -O b \\
                -o "TMP/jeter.bcf" --file-list "TMP/jeter.list"

        bcftools index --force --threads ${task.cpus}  "TMP/jeter.bcf"

	# add the contig to change checksum
	echo "${ctg}" >> TMP/jeter.list
	echo "${bed}" >> TMP/jeter.list
	MD5=`cat TMP/jeter.list  | sha1sum | cut -d ' ' -f1`

        mv TMP/jeter.bcf "concat.\${MD5}.${ctg}.bcf"
        mv TMP/jeter.bcf.csi "concat.\${MD5}.${ctg}.bcf.csi"
"""
}


process LEVEL2 {
tag "${ctg}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
label "process_short"
afterScript "rm -rf TMP"
input:
        tuple val(ctg),path("VCF/*")
output:
	tuple val(ctg),path("concat.*"),emit:output
script:
	def args = "-a --rm-dups all"
	def prefix = "concat"+(ctg.equals(ALL_CONTIGS)?"":"."+ctg)
"""
        mkdir -p TMP
	find VCF -type l \\( -name "*.bcf" -o -name "*.vcf.gz" \\)  | sort > TMP/jeter.list

        bcftools concat --threads ${task.cpus} \\
		${args} \\
                -O b \\
                -o "TMP/jeter.bcf" --file-list "TMP/jeter.list"

        bcftools index --force --threads ${task.cpus}  "TMP/jeter.bcf"


        mv TMP/jeter.bcf "${prefix}.bcf"
        mv TMP/jeter.bcf.csi "${prefix}.bcf.csi"
"""
}
