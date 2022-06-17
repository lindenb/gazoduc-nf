/*

Copyright (c) 2022 Pierre Lindenbaum

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
include {getKeyValue;getModules} from '../../modules/utils/functions.nf'

workflow BCFTOOLS_CONCAT_01 {
	take:
		meta /* params */
		vcfs /* a FILE containing the path to the indexed VCF */
		bed  /* path/to/bed file or empty string */
	main:
		vers_ch = concat_version(meta, vcfs, bed)
		each_list_ch = sqrt_files(meta,vcfs).splitCsv()
		concat0_ch = concat0(meta,each_list_ch,bed)
		concat1_ch = concat1(meta,concat0_ch.vcf.collect())
	emit:
		vcf = concat1_ch.vcf
		version = vers_ch.version
	}


process concat0 {
tag "${file(vcfs).name}"
input:
	val(meta)
	val(vcfs)
	val(bed)
output:
	path("concat.0.bcf"),emit:vcf	
	path("concat.0.bcf.csi"),emit:csi
script:
	"""
	hostname 1>&2
	module load ${getModules("bcftools")}

	bcftools concat --threads ${task.cpus} \
		${bed.isEmpty()?"":"--regions-file \"${bed}\""} |\
		--no-version --allow-overlaps --remove-duplicates \
		-O b -o "concat.0.bcf" --file-list "${vcfs}"

	bcftools index --threads ${task.cpus}  "concat.0.bcf"
	"""
	}

process concat_version {
	executor "local"
	input:
		meta
		vcfs
		bed
	output:
		path("version.xml"),emit:version
	script:
	"""
	cat << EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">concat vcf(s) using bcftools</entry>
		<entry key="vcfs">${vcfs}</entry>
		<entry key="bed">${bed}</entry>
	</properties>
	EOF
	"""
	}

process concat1 {
tag "N=${L.size()}"
input:
	val(meta)
	val(L)
output:
	path("${prefix}concat${suffix}"),emit:vcf	
script:
	prefix = getKeyValue(meta,"prefix","")
	sufffix = getKeyValue(meta,"suffix",".bcf")

	if(L.size()==1)
	"""
	hostname 1>&2
	module load ${getModules("bcftools")}
	
	bcftools sort -T . \
		-O "${suffix.contains("b")?"b":"z"} \
		-o "${prefix}concat${suffix.contains("b")?".bcf":".vcf.gz"}" \
		"${L[0]}"

	bcftools index --threads ${task.cpus} \
		${prefix}concat${suffix.contains("b")?"":"--tbi"} \
		"${prefix}concat${suffix.contains("b")?".bcf":".vcf.gz"}"

	"""
	else
	"""
	hostname 1>&2
	module load ${getModules("bcftools")}
	set -o pipefail

cat << EOF > jeter.list
${L.join(" ")}
EOF

	bcftools concat --threads ${task.cpus} \
		--no-version --allow-overlaps --remove-duplicates \
		-O u --file-list jeter.list |\
		bcftools sort -T . \
			-O "${suffix.contains("b")?"b":"z"} \
			-o "${prefix}concat${suffix.contains("b")?".bcf":".vcf.gz"}" 

	bcftools index --threads ${task.cpus} \
		${prefix}concat${suffix.contains("b")?"":"--tbi"} \
		"${prefix}concat${suffix.contains("b")?".bcf":".vcf.gz"}"

	rm jeter.list
	"""
	}

