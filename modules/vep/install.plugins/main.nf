

process VEP_INSTALL_PLUGINS {
tag "${meta.id}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta),val(plugins)
output:
	tuple val(meta),path("vep.plugins.dir"),emit:directory
	path("versions.yml"),emit:versions
script:
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

mkdir -p vep.plugins.dir

# Install VEP plugins
vep_install --AUTO p --NO_UPDATE --PLUGINSDIR vep.plugins.dir --NO_HTSLIB -g '${plugins}' 1>&2

cat << END_VERSIONS > versions.yml
"${task.process}":
	vep: todo
END_VERSIONS
"""
}

