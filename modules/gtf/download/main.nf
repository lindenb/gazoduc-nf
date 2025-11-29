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
include { parseBoolean } from '../../utils/functions.nf'
include { isBlank      } from '../../utils/functions.nf'


process DOWNLOAD_GTF_OR_GFF3 {
tag "${meta.id}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta),path(dict)
output:
	tuple val(meta),path("*.gz"), path("*.gz.tbi"),optional:true,emit:gtf
	path("versions.yml"),emit:versions
script:
	def extension="";
	def url0= task.ext.url?:""
	def enable_missing = parseBoolean(task.ext.enable_missing?:false)
	if(url0.isEmpty()) {
		 extension = (task.ext==null || task.ext.suffix==null?
				(task.process.toString().toLowerCase().endsWith("gtf")?"gtf":
					(task.process.toString().toLowerCase().endsWith("gff3")?"gff3":""))
				:(task.ext.suffix?:"")
				)
		if(extension.isEmpty()) throw new IllegalArgumentException("suffix missing for ${task.process}");


		def ucsc_name = (task.ext.ucsc_name?:meta.ucsc_name).toString()
		if(isBlank(ucsc_name))  {
			if(!enable_missing) throw new IllegalArgumentException("undefined ucsc_name for ${task.process}");
			url = ""
			}
		else
			{
			def release = task.ext.release?:"114"
			def base = task.ext.base?:"https://ftp.ensembl.org/pub"
			if(ucsc_name=="hg38") {
				url = "${base}/release-${release}/${extension}/homo_sapiens/Homo_sapiens.GRCh38.${release}.${extension}.gz"
				}
			else if(ucsc_name=="hg19") {
				url = "${base}/grch37/current/${extension}/homo_sapiens/Homo_sapiens.GRCh37.87.chr.${extension}.gz"
				}
			else
				{
				if(!enable_missing) throw new IllegalArgumentException("url missing for ${task.process}: '${ucsc_name}'.");
				url = ""
				}
			}
		}
	else
		{
		url = url0;
		extension = url.contains(".gtf")?"gtf":"gff3"
		}
	def prefix = task.ext.prefix?:"${dict.baseName}.${extension}"
"""
hostname 1>&2
mkdir -p TMP
set -o pipefail

if ${!(url.isEmpty() && enable_missing==true)}
then


	curl -L -o TMP/gencode.txt.gz "${url}"

	gunzip -c TMP/gencode.txt.gz |\\
			jvarkit bedrenamechr -f "${dict}" --column 1 --convert SKIP |\\
			LC_ALL=C sort -T TMP -t '\t' -k1,1 -k4,4n |\\
			bgzip > "${prefix}.gz"

	tabix -f -p gff  "${prefix}.gz"

fi

cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(tabix version | awk '(NR==1) {print \$NF;}')"
	URL: "${url}"
END_VERSIONS
"""

stub:
   def prefix = "${dict.baseName}"
   def	 extension = (task.ext==null || task.ext.suffix==null?
				(task.process.toString().toLowerCase().endsWith("gtf")?"gtf":
					(task.process.toString().toLowerCase().endsWith("gff3")?"gff3":""))
				:(task.ext.suffix?:"")
				)
   if(extension.isEmpty()) throw new IllegalArgumentException("suffix missing for ${task.process}");

"""
touch versions.yml "${prefix}.${suffix}.gz" "${prefix}.${suffix}.gz.tbi"
"""
}

