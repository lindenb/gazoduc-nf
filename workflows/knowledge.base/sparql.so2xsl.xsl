<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet
 version="1.0"
 xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
 xmlns:q="http://www.w3.org/2005/sparql-results#"
 xmlns="http://www.w3.org/2000/svg"
 xmlns:axsl="http://www.w3.org/1999/XSL/Transform2"
 >
<xsl:namespace-alias stylesheet-prefix="axsl" result-prefix="xsl"/>

<xsl:output method="xml" encoding="UTF-8" indent="yes"/>

<xsl:template match="/">
<axsl:stylesheet version="1.0">
<xsl:apply-templates select="q:sparql"/>
</axsl:stylesheet>
</xsl:template>

<xsl:template match="q:sparql">

<axsl:template name="biotype2qname">
<axsl:param name="type"/>
<axsl:param name="default"/>
<axsl:choose>
  <xsl:apply-templates select="q:results/q:result"/>
  <axsl:otherwise>
	<axsl:message>undefined biotype <axsl:value-of select="$type"/> for <axsl:value-of select="attributes/attribute[@key='gene_id']"/></axsl:message>
	<axsl:value-of select="$default"/>
  </axsl:otherwise>
</axsl:choose>
</axsl:template>

</xsl:template>

<xsl:template match="q:result">
<axsl:when>
    <xsl:attribute name="test">
    <xsl:text>$type = &apos;</xsl:text>
    <xsl:value-of select="q:binding[@name='label']/q:literal/text()"/>
    <xsl:text>&apos;</xsl:text>
    </xsl:attribute>
    <xsl:text>obo:</xsl:text>
    <xsl:value-of select="substring-after(q:binding[@name='url']/q:uri/text(),'obo/')"/>
</axsl:when>
</xsl:template>

</xsl:stylesheet>