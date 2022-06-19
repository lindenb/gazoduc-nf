process VERSION_TO_HTML {
executor "local"
tag "${xml.name}"
input:
	val(meta)
	path("xml")
output:
	path("${prefix}version.html"),emit:html
script:
	prefix = meta.getOrDefault("prefix","")
"""

cat << __EOF__ > jeter.xsl
<?xml version='1.0' encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl='http://www.w3.org/1999/XSL/Transform' version='1.0'>
<xsl:output method="html"  encoding="UTF-8"/>
<xsl:template match="/">
<html>
<head>
<title>${prefix}</title>
<meta charset="utf-8" />
<style type="text/css">
  dl {
    display: flex;
    flex-flow: row wrap;
    border: solid #333;
    border-width: 3px 1px 0 0;
  }
  dt {
    flex-basis: 20%;
    padding: 2px 4px;
    background: #333;
    text-align: right;
    color: #fff;
  }
  dd {
    flex-basis: 70%;
    flex-grow: 1;
    margin: 0;
    padding: 2px 4px;
    border-bottom: 1px solid #333;
  }

</style>
</head>
<body>
<div>
<xsl:apply-templates select="*"/>
</div>
<hr/>
Pierre Lindenbaum PhD. Institut du Thorax. 44400 Nantes. France.
</body>
</html>
</xsl:template>


<xsl:template match="properties">
<dl>
<xsl:apply-templates select="*"/>
</dl>
</xsl:template>

<xsl:template match="entry">
<dt><xsl:value-of select="@key"/></dt>
<dd><xsl:apply-templates></dd>
</xsl:template>


</xsl:stylesheet>
__EOF__


xsltproc jeter.xsl "${xml}" > "${prefix}version.html"

rm jeter.xsl
"""
}

