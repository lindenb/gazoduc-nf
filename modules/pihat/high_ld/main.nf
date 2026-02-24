/*
Copyright (c) 2026 Pierre Lindenbaum

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
include { verify       } from '../../utils/functions.nf'

/*
 * Download regions with HIGH LD 
 * see : https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)#cite_note-3
 * https://github.com/meyer-lab-cshl/plinkQC/blob/master/plinkQC.Rcheck/plinkQC/extdata/high-LD-regions-hg38-GRCh38.txt
 *
 */
process DOWNLOAD_HIGH_LD {
tag "${meta.id}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta),path(dict)
output:
	tuple val(meta),path("*.bed"),emit:bed
	path("versions.yml"),emit:versions
script:
	def extension="";
	def url= task.ext.url?:""
	if(isBlank(url)) {
		if(meta.ucsc_name=="hg38") {
			url = "https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/refs/heads/master/plinkQC.Rcheck/plinkQC/extdata/high-LD-regions-hg38-GRCh38.txt"
			}
		else if(meta.ucsc_name=="hg19") {
			url = "https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/refs/heads/master/plinkQC.Rcheck/plinkQC/extdata/high-LD-regions-hg19-GRCh37.bed"
			}
		}
	def prefix = task.ext.prefix?:"${dict.baseName}.highld"
"""
hostname 1>&2
mkdir -p TMP
set -o pipefail

touch "${prefix}.bed"

if ${!isBlank(url)}
then

	curl -L "${url}" |\\
		awk '{printf("%s\t%s\t%s\\n",\$1,\$2,\$3);}' |\\
		jvarkit bedrenamechr -f "${dict}" --column 1 --convert SKIP |\\
		LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n > "${prefix}.bed"


fi

cat << END_VERSIONS > versions.yml
"${task.process}":
	URL: "${url}"
END_VERSIONS
"""

stub:
	def prefix = task.ext.prefix?:"${dict.baseName}.highld"
"""
touch versions.yml "${prefix}.bed"
"""
}

