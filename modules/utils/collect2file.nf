process collect2file {
executor "local"
tag "N=${L.size()}"
input:
	val(meta)
	val(L)
output:
	path("concat.list"),emit:output
script:
"""
hostaname 1>&2
set -o pipefail

cat << EOF | awk -F '/' '{printf("%s\t%s\\n",\$NF,\$0);}' | sort -t '\t' -T. -k1,1 -k2,2 | cut -f 2 | uniq > concat.list
${L.join("\n")}
EOF
"""
}
