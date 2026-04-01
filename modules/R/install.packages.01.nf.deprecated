include {moduleLoad;getKeyValue;assertNotEmpty} from '../utils/functions.nf'

process R_INSTALL_PACKAGES_01 {
input:
	val(meta)
output:
        path("LIB"),emit:lib
	path("version.xml"),emit:version
script:
	def r_version= getKeyValue(meta,"R_module","r/3.6.3")
	def r_repos =  getKeyValue(meta,"R_repos","http://cran.r-project.org")
	def r_packages =  assertNotEmpty(getKeyValue(meta,"R_packages",""),"R_packages must be declared")
"""
hostname 1>&2
${moduleLoad(r_version)}

mkdir LIB
cat << EOF | R --vanilla
install.packages(${r_packages},repos="${r_repos}",lib="\${PWD}/LIB")
EOF

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">install package(s) for R</entry>
        <entry key="R.repos">${r_repos}</entry>
        <entry key="R.packages">${r_packages}</entry>
	<entry key="R.version">\$(R --version | head -n1)</entry>
</properties>
EOF
"""
stub:
"""
mkdir LIB
echo "<properties/>" > version.xml
"""
}

