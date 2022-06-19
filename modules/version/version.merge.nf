process MERGE_VERSION {
tag "N=${L.size()}"
executor "local"
input:
	val(meta)
	val(name)
	val(description)
	val(L)
output:
	path("version.xml"),emit:version
script:
	prefix = meta.getOrDefault("prefix","")
"""
cat << EOF > jeter.xml
<properties id="${name}">
	<entry key="name">${name}</entry>
	<entry key="description">${description}</entry>
	<entry key="date">\$(date)</entry>
	<entry key="steps">
EOF

for X in ${L.join(" ")}
do
	xmllint --format "\${X}" | tail -n+2 >> jeter.xml
done

cat << EOF >> jeter.xml
	</entry>
</properties>
EOF
	xmllint --format jeter.xml > "${prefix}version.xml"
rm jeter.xml
"""
}
