def ALL_CONTIGS="ALL"

List makeSQRT(def L1) {
        def key = L1.get(0);
        def L = L1.get(1);
        int n = (int)Math.ceil(Math.sqrt(L.size()));
        if(n<25) n=25;
        def returnList = [];
        def currList = [];
        int i=0;
	//System.err.println("L.size "+L.size()+ " "+n+" "+L);
        for(;;) {
                if(i<L.size()) currList.add(L.get(i));
                if(i==L.size() || currList.size()==n) {
                        if(!currList.isEmpty()) returnList.add([key,currList]);
                        if(i==L.size()) break;
                        currList=[];
                        }
                i++;
                }
	//System.err.println("return "+returnList);
        return returnList;
        }

workflow BCFTOOLS_CONCAT_CONTIGS {
take:
	vcfs // tuple [contig,vcf,idx]
	bed //optional bed file
main:
	ch2 = vcfs.
		map{[it[0],[it[1],it[2]]]}.
		groupTuple().
		map{[it[0],it[1].sort{t->t[0].name}]}. // prevent cache invalidation due to order
		flatMap{makeSQRT(it)}.
		map{[it[0], it[1].flatten()]}

	ch3 = LEVEL1(ch2,bed)
	ch4 = LEVEL2(ch3.output.groupTuple().map{[it[0],it[1].flatten().sort()]})

emit:	
	output = ch4.output

}


process LEVEL1 {
tag "${ctg} ${bed.name}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
label "process_quick"
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
