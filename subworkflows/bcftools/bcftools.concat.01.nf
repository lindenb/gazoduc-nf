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
include {isBlank;moduleLoad;getKeyValue;getModules;getClassTaxonomy} from '../../modules/utils/functions.nf'
include {SQRT_FILE} from '../../modules/utils/sqrt.nf'


String regionsArgs(def row) {
	if(!isBlank(row.bed)) return "--regions-file \"${row.bed}\" ";
	if(!isBlank(row.contig)) return "--regions \"${row.contig}\" ";
	if(!isBlank(row.interval)) return "--regions \"${row.interval}\" ";
	return "";
	}

workflow BCFTOOLS_CONCAT_01 {
take:
	meta
	row
main:
	def optR = regionsArgs(row)
	d1_ch = SQRT_FILE(meta, row.vcfs)

	d2_ch = d1_ch.clusters.splitText().map{it.trim()}
	d3_ch = BCFTOOL_CONCAT_01(meta, optR, d2_ch )
	d4_ch = BCFTOOL_CONCAT_02(meta, optR, d3_ch.vcf.collect() )
emit:
	vcf = d4_ch.vcf
}


process BCFTOOL_CONCAT_01 {
tag "${file(vcfs).name}"
cpus 1
input:
        val(meta)
	val(optR)
        val(vcfs)
output:
        path("concat.0.bcf"),emit:vcf
	path("version.xml"),emit:version
script:
"""

	hostname 1>&2
	${moduleLoad("bcftools")}

	bcftools concat --threads ${task.cpus} ${optR} \
		--no-version --allow-overlaps --remove-duplicates \
		-O b -o "concat.0.bcf" --file-list "${vcfs}"

	bcftools index --threads ${task.cpus}  "concat.0.bcf"

	touch version.xml
"""
}

process BCFTOOL_CONCAT_02 {
input:
        val(meta)
	val(optR)
        val(vcfs)
output:
        path("concat.1.bcf"),emit:vcf
	path("version.xml"),emit:version
script:
"""

	hostname 1>&2
	${moduleLoad("bcftools")}

cat << EOF > jeter.list
${vcfs.join("\n")}
EOF

	bcftools concat --threads ${task.cpus} ${optR} \
		--no-version --allow-overlaps --remove-duplicates \
		-O b -o "concat.1.bcf" --file-list jeter.list

	bcftools index --threads ${task.cpus}  "concat.1.bcf"

	touch version.xml

"""
}



workflow __BCFTOOLS_CONCAT_01 {
	take:
		meta

		/* row.vcfs a FILE containing the path to the indexed VCF */
		/* row.contig or empty string */
		/* row.interval or empty string */
		/* row.bed path/to/bed file or empty string */	
		row
	main:
		//if(!row.containsKey("vcfs")) throw new IllegalArgumentException("undefined row.vcfs in "+row);
		version_ch = CONCAT_VERSION(meta,row)

		each_list_ch = SQRT_FILE(meta, row.vcfs)
		concat0_ch = CONCAT0(meta, each_list_ch.clusters.splitText().map{it.trim()}, row)
		concat1_ch = CONCAT1(meta,concat0_ch.vcf.collect())
	emit:
		vcf = concat1_ch.vcf
		index = concat1_ch.index
		version = version_ch.version

	}


process CONCAT0 {
tag "${vcfs}"
input:
	val(meta)
	val(vcfs)
	val(row)
output:
	path("concat.0.bcf"),emit:vcf	
	path("concat.0.bcf.csi"),emit:csi
script:
	def optR = regionsArgs(row)
	"""
	hostname 1>&2
	${moduleLoad("bcftools")}

	bcftools concat --threads ${task.cpus} ${optR} \
		--no-version --allow-overlaps --remove-duplicates \
		-O b -o "concat.0.bcf" --file-list "${vcfs}"

	bcftools index --threads ${task.cpus}  "concat.0.bcf"
	"""
	}

process CONCAT_VERSION {
	executor "local"
	input:
		val(meta)
		val(row)
	output:
		path("version.xml"),emit:version
	script:
	"""
	hostname 1>&2
	${moduleLoad("bcftools")}

	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">concat vcf(s) using bcftools</entry>
		<entry key="vcfs">${row.vcfs}</entry>
		<entry key="count(vcfs)">\$(wc -l < ${row.vcfs})</entry>
		<entry key="bed">${row.bed?:""}</entry>
		<entry key="contig">${row.contig?:""}</entry>
		<entry key="interval">${row.interval?:""}</entry>
		<entry key="bcftools">\$( bcftools --version-only)</entry>
	</properties>
	EOF
	"""
	}

process CONCAT1 {
tag "N=${L.size()}"
input:
	val(meta)
	val(L)
output:
	path("${prefix}concat${suffix}"),emit:vcf	
	path("${prefix}concat${suffix}${suffix.contains("b")?".csi":".tbi"}"),emit:index
script:
	prefix = getKeyValue(meta,"prefix","")
	suffix = getKeyValue(meta,"suffix",".bcf")

	if(L.size()==1)
	"""
	hostname 1>&2
	${moduleLoad("bcftools")}
	
	bcftools sort -T . \
		-O "${suffix.contains("b")?"b":"z"}" \
		-o "${prefix}concat${suffix.contains("b")?".bcf":".vcf.gz"}" \
		"${L[0]}"

	bcftools index --threads ${task.cpus} \
		${suffix.contains("b")?"":"--tbi"} \
		"${prefix}concat${suffix.contains("b")?".bcf":".vcf.gz"}"

	"""
	else
	"""
	hostname 1>&2
	${moduleLoad("bcftools")}
	set -o pipefail

cat << EOF > jeter.list
${L.join("\n")}
EOF

	bcftools concat --threads ${task.cpus} \
		--no-version --allow-overlaps --remove-duplicates \
		-O u --file-list jeter.list |\
		bcftools sort -T . \
			-O "${suffix.contains("b")?"b":"z"}" \
			-o "${prefix}concat${suffix.contains("b")?".bcf":".vcf.gz"}" 

	bcftools index --threads ${task.cpus} \
		 ${suffix.contains("b")?"":"--tbi"} \
		"${prefix}concat${suffix.contains("b")?".bcf":".vcf.gz"}"

	rm jeter.list
	"""
	}

