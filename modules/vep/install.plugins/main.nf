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
process VEP_INSTALL_PLUGINS {
tag "${meta.id?:""}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	val(meta)
output:
	tuple val(meta),path("*.plugins.dir"),emit:directory
	path("versions.yml"),emit:versions
script:
	def plugins=task.ext.plugins?:""
	def prefix = task.ext.prefix?:"vep"
"""
set -x
set -o pipefail
VEPPATH=`which vep`
VEPDIR=`dirname \${VEPPATH}`
ls "\${VEPDIR}/../lib/perl5/5.32/core_perl/Config_heavy.pl"
sed -ri 's/([^\\@])@univ-nantes/\\1\\\\@univ-nantes/g' "\${VEPDIR}/../lib/perl5/5.32/core_perl/Config_heavy.pl"

# test user
echo "\${USER}" 1>&2
test ! -z "\${USER}"

mkdir -p "${prefix}.plugins.dir"

if ${!plugins.isEmpty()}
then
	# Install VEP plugins
	vep_install --AUTO p --NO_UPDATE --PLUGINSDIR "${prefix}.plugins.dir" --NO_HTSLIB -g '${plugins}' 1>&2
fi

cat << END_VERSIONS > versions.yml
"${task.process}":
	vep: todo
END_VERSIONS
"""

stub:
"""
mkdir vep.plugins.dir
touch versions.yml
"""
}
