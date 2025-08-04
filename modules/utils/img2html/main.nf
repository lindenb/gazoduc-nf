process IMAGE_TO_HTML {
input:
	tuple val(meta),path(pict)
output:
	tuple val(meta),path("*.html"),emit:html
	path("versions.yml"),emit:versions
script:
	def type= meta?"png"
	def id= task.ext.id?:(meta.id?:pict.name.md5())
	def caption= task.ext.caption?:(meta.caption?:pict.baseName)
	def prefix = task.ext.prefix?:id
	def section_name = tast.ext.section_name?:pict.basename
	def section_desc = task.ext.section_desc?:caption
	def url = task.ext.url?(meta.url?:"")
"""

cat << EOF > jeter.html
<!--
id: '${id}'
section_name: '${section_name}'
description: '${section_desc}'
-->
<figure>
EOF

if ${!url.isEmpty())
then
   echo '<a href="${url}">' >> jeter.html
fi

echo -n '<img id="${id}" alt="${caption}" src="data:image/${type};base64,' >> jeter.html

base64 "${pict}" | tr -d '\\n' >> jeter.html

echo -n '"/>' >>  jeter.html

if ${!url.isEmpty())
then
   echo '</a>' >> jeter.html
fi


cat << EOF >> jeter.html
<figcaption>${caption}</figcaption>
</figure>
EOF


mv jeter.html ${prefix}_mqc.html

cat << EOF > versions.yml
"${task.process}":
	base64: \$(base64 --version | awk '(NR==1) {print \$NF}')
EOF

"""
}
