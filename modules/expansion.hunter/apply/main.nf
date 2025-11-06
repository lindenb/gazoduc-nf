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
process XHUNTER_APPLY {
	label "process_single"
	tag "${meta.id}"
	afterScript "rm -rf TMP"
    conda "${moduleDir}/../../../conda/xhunter.yml"
	input:
		tuple val(meta1),path(fasta)
		tuple val(meta2),path(fai)
		tuple val(meta3),path(catalog)
		tuple val(meta ),path(bam),path(bai)
	output:
		tuple val(meta),path("*.vcf.gz"),path("*.tbi"),emit:vcf
		tuple val(meta),path("*.json"),optional:true,emit:json
		tuple val(meta),path("*.bam"),path("*.bai"),optional:true,emit:bam
		path("versions.yml"),emit:versions
	script:
		def keep_bams = (task.ext.keep_bams?:false).toBoolean()
		def keep_json = (task.ext.keep_json?:false).toBoolean()
		def mode      = (task.ext.mode?:"streaming")
		def prefix    = task.ext.prefix?:"${meta.id}"
		def sex       = meta.sex?:"undefined"
		def args1     = task.ext.args1?:""
	
	//log.warn("keep_bams = ${keep_bams} ${task.ext.keep_bams} class2=${task.ext.keep_bams.class} class=${keep_bams.class}")
	"""
	hostname 1>&2
	mkdir -p TMP
	export TMPDIR=\${PWD}/TMP
	
	ExpansionHunter \\
		${args1} \\
		--reads "${bam}" \\
		--reference "${fasta}" \\
		--variant-catalog "${catalog}" \\
		--threads ${task.cpus} \\
		--output-prefix "TMP/jeter" \\
		--analysis-mode "${mode}" \\
		--sex '${sex}' \\
		--log-level info 1>&2

	# add missing INFO in header
	awk '
	/^##/    {
		if(\$0 ~ /INFO=<ID=END,/) info_END = 1;
		else if(\$0 ~ /INFO=<ID=REF,/) info_REF = 1;
		else if(\$0 ~ /INFO=<ID=PL,/) info_PL = 1;
		else if(\$0 ~ /INFO=<ID=RU,/) info_RU = 1;
		else if(\$0 ~ /INFO=<ID=REPID,/) info_REPID = 1;
		else if(\$0 ~ /FORMAT=<ID=SO,/) fmt_SO = 1;
		else if(\$0 ~ /FORMAT=<ID=REPCN,/) fmt_REPCN = 1;
		else if(\$0 ~ /FORMAT=<ID=REPCI,/) fmt_REPCI = 1;
		else if(\$0 ~ /FORMAT=<ID=ADSP,/) fmt_ADSP = 1;
		else if(\$0 ~ /FORMAT=<ID=ADFL,/) fmt_ADFL = 1;
		else if(\$0 ~ /FORMAT=<ID=ADIR,/) fmt_ADIR = 1;
		}

	/^#CHROM/ {
		if(info_END!=1)    printf("##INFO=<ID=END,Number=1,Type=Integer,Description=\\"END variant\\">\\n");
		if(info_REF!=1)    printf("##INFO=<ID=REF,Number=1,Type=Integer,Description=\\"missing desc\\">\\n");
		if(info_PL!=1)     printf("##INFO=<ID=RL,Number=1,Type=Integer,Description=\\"missing desc\\">\\n");
		if(info_RU!=1)     printf("##INFO=<ID=RU,Number=1,Type=String,Description=\\"missing desc\\">\\n");
		if(info_REPID!=1)  printf("##INFO=<ID=REPID,Number=1,Type=String,Description=\\"missing desc\\">\\n");


		if(fmt_SO!=1)    printf("##FORMAT=<ID=SO,Number=1,Type=String,Description=\\"missing desc\\">\\n");
		if(fmt_REPCN!=1) printf("##FORMAT=<ID=REPCN,Number=1,Type=String,Description=\\"Number of repeat units spanned by the allele\\">\\n");
		if(fmt_REPCI!=1) printf("##FORMAT=<ID=REPCI,Number=1,Type=String,Description=\\"missing desc\\">\\n");
		if(fmt_ADSP!=1)  printf("##FORMAT=<ID=ADSP,Number=1,Type=String,Description=\\"missing desc\\">\\n");
		if(fmt_ADFL!=1)  printf("##FORMAT=<ID=ADFL,Number=1,Type=String,Description=\\"missing desc\\">\\n");
		if(fmt_ADIR!=1)  printf("##FORMAT=<ID=ADIR,Number=1,Type=String,Description=\\"missing desc\\">\\n");
		}
		{print}' TMP/jeter.vcf > TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter.vcf
	
	# rename sample
	bcftools query -l TMP/jeter.vcf | awk '{printf("%s\t${meta.id}\\n",\$1);}' > TMP/rename.txt

	bcftools reheader --fai "${fai}" --samples TMP/rename.txt TMP/jeter.vcf |\\
		awk '/^#CHROM/ {printf("##samples.sex=${meta.id}=${meta.sex}\\n",S);} {print;}' > TMP/jeter2.vcf
	
	mv TMP/jeter2.vcf TMP/jeter.vcf
	
	#sort and index
	bcftools annotate --set-id '%INFO/REPID' -O u TMP/jeter.vcf |\\
		bcftools sort --max-mem '${task.memory.giga}G' -T TMP/sort -O z -o "TMP/${prefix}.vcf.gz"
	bcftools index   --threads ${task.cpus} -f -t "TMP/${prefix}.vcf.gz"


	if  ${keep_bams}
	then
		samtools sort -T TMP/sort  --threads ${task.cpus} -O BAM --reference "${fasta}" -o "TMP/${prefix}.realigned.bam" "TMP/jeter_realigned.bam"
		samtools index --threads ${task.cpus} "TMP/${prefix}.realigned.bam"
		mv "TMP/${prefix}.realigned.bam" ./
		mv "TMP/${prefix}.realigned.bam.bai" ./
	fi


	mv TMP/${prefix}.vcf.gz ./
	mv TMP/${prefix}.vcf.gz.tbi ./

	if ${keep_json}
	then
		mv TMP/jeter.json "${prefix}.json"
	fi

cat << EOF > versions.yml
"${task.process}":
    xhunter: \$( ExpansionHunter --version 2>&1 | awk '/^Start/ {print \$NF}')
EOF
"""
stub:
"""
touch versions.yml "${meta.id}.json" "${meta.id}.vcf.gz" "${meta.id}.vcf.gz.tbi" "${meta.id}.realigned.bam" "${meta.id}.realigned.bam.bai"
"""
}
