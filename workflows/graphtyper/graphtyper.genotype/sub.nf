String toFilePath (String contig,String start,String end) {
	def dir2 = contig+":"+start+"-"+end;
	dir2 = dir2.md5();
	dir2 = ""+workflow.workDir + "/TMP_VCFS/" + contig + "/" + dir2.substring(0,2) + "/" + dir2.substring(2,4) + "/"+ params.prefix + contig + "_" + ((start as int)-1) +"_" + end + ".bcf";
	return dir2;
	}


String toFailedPath(String contig,String start,String end) {
	return toFilePath(contig,start,end)+".failed";
	}


def touchFailed(String contig,String start,String end) {
	def f1 = new java.io.File(toFailedPath(contig,start,end));
	System.err.println("Create FAILED file "+f1);
	f1.getParentFile().mkdirs();
	f1.createNewFile();
	return [contig,start,end];
	}

def splitIntervals(chrom,start1,end1) {
	try {
	int start = (start1 as int) - 1;
	int end = (end1 as int);
	int length = end-start ;
	int l2 = Math.max(1,(int)(length/10));
	def L = [];
	if(l2>=length) {System.err.println("Skipping ${chrom}:${start1}-${end1}") ; return L;}
	while(start<end) {
		L.add([chrom , String.valueOf(start+1) , String.valueOf(Math.min(start+l2,end))]);
		start+=l2;
		}
	return L;
	} catch(Throwable err) {
		System.err.println("BOUM in splitIntervals '${chrom}:${start1}-${end1}'");
		err.printStackTrace();
		throw err;
		}
	}


workflow DIVIDE_AND_CONQUER {
	take:
		level
		samplesheet
		intervals
	main:
		ok_ch = Channel.empty()

		intervals.branch {
			done : file(toFilePath(it[0],it[1],it[2])).exists() && file(toFilePath(it[0],it[1],it[2])+".csi").exists()
			known_to_fail: file(toFailedPath(it[0],it[1],it[2])).exists()
			todo: true
			}.set{branch1}


		//branch1.done.view{"ALREADY DONE: ${it[0]+":"+it[1]+"-"+it[2]}"}

		
		ok_ch = ok_ch.mix(branch1.done.map{[ it[0], file(toFilePath(it[0],it[1],it[2])) ]})


		sentinel_ch = branch1.todo.map{[ [it[0], it[1], it[2] ], file("NO_VCF")]}
        	call_ch = GRAPHTYPER(level, samplesheet, branch1.todo)

        	call_ch.output.map{[ [ it[0], it[1], it[2] ], it[3] ]}.
			mix(sentinel_ch).
			groupTuple().
        	        branch {
                	        failed : it[1].size()==1
                        	ok: true
	                }.set{result_ch}


		//result_ch.ok.view{"OK: ${it[0][0]+":"+it[0][1]+"-"+it[0][2]}"}

		ok_ch = ok_ch.mix( result_ch.ok.map{[it[0][0], file(toFilePath(it[0][0],it[0][1],it[0][2])) ]} )

		failed1_ch = result_ch.failed.
			map{[it[0][0], it[0][1], it[0][2]]} 

		failed1_ch.
			map{touchFailed(it[0],it[1],it[2])}.
			view{"FAILED: ${it[0]+":"+it[1]+"-"+it[2]}"}

				

		failed_ch= failed1_ch.
			mix(branch1.known_to_fail).
			flatMap{splitIntervals(it[0],it[1],it[2])}

	emit:
		ok = ok_ch
		failed= failed_ch
	}

boolean isClusterError(task) {
	if(task==null) return false;
	if(task.previousException == null) return false;
	String s = task.previousException.message;
	if(s.startsWith("Error submitting process ") && s.endsWith("for execution")) {
		return true;
		}
	return false;
	}

process GRAPHTYPER {
label "process_single"
tag "(level ${level}) ${contig}:${start}-${end}"
conda "${moduleDir}/../../../conda/graphtyper.yml"
afterScript "rm -rf TMP TMP2"
errorStrategy  = {isClusterError(task)?"terminate":"ignore"}
input:
	val(level)
	path(samplesheet)
	tuple val(contig),val(start),val(end)
output:
	tuple val(contig), val(start),val(end), path("*.flag"),emit:output
script:
	def outputname = toFilePath(contig,start,end)
	def copy_bams = (level > 3)
	


if(file(outputname).exists() && file(outputname+".csi").exists())
"""
echo "${outputname}" > "${contig}_${start}_${end}.flag"
"""
else
"""
module  load jvarkit
mkdir -p TMP
mkdir -p TMP2
export TMPDIR=\${PWD}/TMP
set -x


cut -f 3 "${samplesheet}" > TMP/bams.list
awk -F '\t' '{print ((\$5*.1.0)/(\$2*1.0));}' "${samplesheet}"  > TMP/cov_by_readlen.txt


if ${copy_bams} ; then

	mkdir -p TMP/BAMS
	# don't forget to erase previous original list....
	rm TMP/bams.list
	i=1
	cut -f 3 "${samplesheet}" | while read SAM
	do
		echo "\$i : \${SAM}" 1>&2

		# extract reads in regions
		samtools view --uncompressed -O BAM  -T "${params.fasta}" "\${SAM}" "${contig}:${start}-${end}" |\\
			java -jar \${HOME}/jvarkit.jar samrmdupnames \\
				-R "${params.fasta}" \\
				--validation-stringency SILENT \\
				--samoutputformat BAM \\
				--bamcompression 0 |\\
			samtools view --reference "${params.fasta}" \\
				 -O BAM -o TMP/BAMS/tmp.jeter.bam
		
		samtools index -@ "${task.cpus}" TMP/BAMS/tmp.jeter.bam

		mv -v "TMP/BAMS/tmp.jeter.bam" "TMP/BAMS/tmp.\${i}.bam"
		mv -v "TMP/BAMS/tmp.jeter.bam.bai" "TMP/BAMS/tmp.\${i}.bam.bai"

		echo "TMP/BAMS/tmp.\${i}.bam" >> TMP/bams.list

		i=\$((i+1))

	done

	test -s TMP/bams.list
fi



graphtyper genotype \\
	"${params.fasta}" \\
	--output=TMP2 \\
	--force_no_copy_reference \\
	--force_use_input_ref_for_cram_reading \\
	--sams=TMP/bams.list \\
	--region="${contig}:${start}-${end}" \\
	--threads=${task.cpus}

rm -rf "TMP2/input_sites"

find TMP2 -type f -name "*.vcf.gz"  | while read F
do
	bcftools query -l "\${F}" | LC_ALL=C sort -T .  > TMP/ordered.samples.txt
	test -s TMP/ordered.samples.txt
	
	bcftools view  --threads ${task.cpus} -O b -o "\${F}.bcf"  --samples-file TMP/ordered.samples.txt "\${F}"
	bcftools index --threads ${task.cpus} "\${F}.bcf"
	echo "\${F}.bcf" >> TMP/merge.bcf.list
	rm "\${F}" "\${F}.tbi"
done


bcftools concat -a --file-list TMP/merge.bcf.list --threads ${task.cpus} -O b9 --write-index -o TMP/merged.bcf

mkdir -p  `dirname "${outputname}"`
mv -v TMP/merged.bcf "${outputname}"
mv -v TMP/merged.bcf.csi "${outputname}.csi"

echo "${outputname}" > "${contig}_${start}_${end}.flag"
"""
}
