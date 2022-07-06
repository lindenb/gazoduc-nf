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

process VERSION_TO_HTML {
executor "local"
tag "${xml.name}"
input:
	val(meta)
	path(xml)
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
${params.rsrc.author}
<br/>
<img src="data:image/jpg;base64,/9j/4AAQSkZJRgABAQAAAQABAAD/2wCEAAkGBxMQEhUQExMWFhUVFxUVEBYXGBUZFhYVFRUdFxcXGxoZHSggGBomGxgXITEhMSorLi4uGB8zODMtNygtLisBCgoKDg0OGxAQGy8lICYtLS4tLS0vLys1LS0tLy4tLS0tLSstLS0tLS0tLS0tLS0tLS0tLS8tKy0tLSstLS0tLf/AABEIAMsA+QMBIgACEQEDEQH/xAAcAAEAAgMBAQEAAAAAAAAAAAAABgcDBAUCAQj/xABJEAABAwIEAgUHCQUGBQUAAAABAAIDBBEFEiExBkEHEyJRcRQyYYGRobEjNDVCUnJzsrNigpLB0SUzQ3TC4RUkU6LSFjZVg8P/xAAbAQEAAwEBAQEAAAAAAAAAAAAAAgMEAQUGB//EADURAAEDAQYDBwQBBAMBAAAAAAEAAhEDBBIhMUFRYXHBE4GRobHR8AUiMjPhI0LC8WJyshT/2gAMAwEAAhEDEQA/ALxRERERERERERERERERERERERERFAuO+L6ihmbFC2Mgsa8l7XE3LnDk4adn3qynTdUddbmq6tVtNt52SnqLhcGYs+to4qmQND35w4NuG9iRzLgEki+W+67qg5paSDopgyJRERcXUREREREREREREREREREREREREREREREREREREREREREREVacedIHVF1NRuBeLtlmGoYdi1ne/vOw9J2+9JXGRjvQ07rPItUSA6sB/wANp5OI3PIek6RHgvg9+IPubsgYbSPG5P2GftW58vYDuoUGtb2tXLT56BZKtYl1xman3RBUvkonl73PIneAXOLiBkYbXOu5J9anS08OoI6eNsUTAxjRZrR7z6SdydytxZKrw95cNVpY260BFUnS785Z+C3871baqTpd+cs/Bb+d602H9vcVmt36TzHqpf0W/RkHjN+u9dXiTGG0cDpyLkENY29sznHQXsbaXPgCuV0W/RkHjN+u9RnpWxXNIymadGDO/wC+4aDxDfzrLaXXXvPE+q9b6XZf/orMYcszyHuYHepXwvxbFXExhrmSNbmc02c2wIBIcN9SNwCpOoL0WYV1dO6ocO1M45fw2Ege12Y+Fl2uKuI2UEYeW5nuNo2XtfvJNjYAejcgc1U1xuS5WWqztNrNGzjWAJ11xOkznoJlSBFw+G+Ioq5hfHdrmW6xjtxmvbXZwNjY+jku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"/>
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
<dd><xsl:apply-templates/></dd>
</xsl:template>


<xsl:template match="pre">
<pre><xsl:apply-templates/></pre>
</xsl:template>

<xsl:template match="b">
<b><xsl:apply-templates/></b>
</xsl:template>

<xsl:template match="i">
<i><xsl:apply-templates/></i>
</xsl:template>

<xsl:template match="code">
<code><xsl:apply-templates select="*"/></code>
</xsl:template>

<xsl:template match="a">
<xsl:element name="a">
<xsl:attribute name="href">
<xsl:choose>
	<xsl:when test="@href"><xsl:value-of select="@href"/></xsl:when>
	<xsl:otherwise><xsl:value-of select="."/></xsl:otherwise>
</xsl:choose>
</xsl:attribute>
<xsl:apply-templates/>
</xsl:element>
</xsl:template>

<xsl:template match="img">
<xsl:element name="img">
<xsl:attribute name="title">
	<xsl:value-of select="."/>
</xsl:attribute>
<xsl:attribute name="src">
<xsl:choose>
	<xsl:when test="@src"><xsl:value-of select="@src"/></xsl:when>
	<xsl:otherwise><xsl:value-of select="."/></xsl:otherwise>
</xsl:choose>
</xsl:attribute>
</xsl:element>
</xsl:template>

<xsl:template match="text()">
<xsl:value-of select="."/>
</xsl:template>

</xsl:stylesheet>
__EOF__


xsltproc jeter.xsl "${xml}" > "${prefix}version.html"

rm jeter.xsl
"""

stub:
"""
touch "${prefix}version.html"
"""
}

