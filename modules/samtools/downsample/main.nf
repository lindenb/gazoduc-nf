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

process DOWNSAMPLE {
   tag "${meta.id} / ${bam.name} DP=${depth}"
   label "process_single"
   afterScript "rm -rf TMP"
   conda "${moduleDir}/../../../conda/bioinfo.01.yml"
   input:
		tuple val(meta1),path(fasta)
		tuple val(meta2),path(fai)
		tuple val(meta3),path(dict)
		tuple val(meta4),path(bed)
	    tuple val(meta),path(bam),path(mosdepth_summary),val(depth)
   output:
		tuple val(meta),path("*.bam"),path("*.bai"),emit:bam
		path("versions.yml"),emit:versions
   script:
		def dps = String.format("%03d",depth)
		def mapq = task.ext.mapq?:1
   """
   hostname 1>&2
  
   mkdir -p TMP

   cut -f 1  "${bed}" | uniq | sort | uniq | while read C
   do
	FRAC=`awk -F '\t' -vC=\$C '(\$1==sprintf("%s_region",C)) {D=(\$4*1.0);if(D==0) D=1E-6; F=(${depth}*1.0)/D ; if(F>0 && F<1.0) print F;}' '${mosdepth_summary}'  `

	echo "${sample} \${C} FRAC=\${FRAC}" >&2 

	awk -F '\t' -vC=\$C '(\$1==C)' "${bed}" > TMP/jeter.bed


	if [ "\${FRAC}" != "" ] ; then

   		samtools view -M -L TMP/jeter.bed --threads ${task.cpus} -q ${mapq} -F 3844 -u --reference "${fasta}" "${bam}"  |\
			samtools view -s \${FRAC} -O BAM -o "TMP/_chunck.\${C}.bam"

	else

   		samtools view -M -L TMP/jeter.bed --threads ${task.cpus} -q ${mapq} -F 3844 -u --reference "${fasta}" \
			-O BAM -o "TMP/_chunck.\${C}.bam" "${bam}"
      
	fi
   done

   samtools merge --threads ${task.cpus} --reference ${fasta} TMP/merged.bam TMP/_chunck.*.bam

   gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" AddOrReplaceReadGroups \
	     I=TMP/merged.bam \
	     O="${sample}.DP${depth}.bam" \
	     RGID=${sample}.DP${dps} \
	     RGLB=${sample}.DP${dps} \
	     RGSM=${sample}.DP${dps} \
	     RGPL=ILLUMINA \
	     RGPU=unit1

   samtools index -@ ${task.cpus} "${sample}.DP${depth}.bam"

   ${mosdepth.toRealPath()} --no-per-base -t ${task.cpus} --by ${bed} --fasta ${reference} \
		--mapq ${params.mapq} TMP/${sample}.after.DP${depth} ${sample}.DP${depth}.bam

   mv "TMP/${sample}.after.DP${depth}.mosdepth.summary.txt" ./

cat << EOF > versions.yml
"${task.process}":
    samtools:\$(samtools  --version | head -n 1| cut -d ' ' -f2)
EOF
   """
   
   stub:
   """
   touch versions.yml ${prefix}.bam ${prefix}.bam.bai
   """
   }
