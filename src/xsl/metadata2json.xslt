<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet
        version="1.0"
        xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
 >

<xsl:output method="text" encoding="UTF-8" />

<xsl:template match="/">
<xsl:apply-templates select="schema"/>
</xsl:template>

<xsl:template match="schema">
{
  "title" : 
  <xsl:choose>
	<xsl:when test="title">
		<xsl:apply-templates select="title" mode="quote"/>
	</xsl:when>
	<xsl:otherwise>
		<xsl:text>"undefined workflow"</xsl:text>
	</xsl:otherwise>
  <xsl:choose>,
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "description":
  <xsl:choose>
	<xsl:when test="description">
		<xsl:apply-templates select="description" mode="quote"/>
	</xsl:when>
	<xsl:otherwise>
		<xsl:text>"no description"</xsl:text>
	</xsl:otherwise>
  <xsl:choose>,
  "type": "object"
  
}
</xsl:template>

<xsl:template match="*" mode="quote">
<xsl:text>"</xsl:text>
<xsl:call-template name="escape_json">
  <xsl:with-param name="s" select="./text()"/>
</xsl:call-template>
<xsl:text>"</xsl:text>
</xsl:template>

<xsl:template name="quote">
<xsl:param name="s"/>
<xsl:value-of select="$s"/>
</xsl:template>


{
  "title": "Graphtyper pipeline parameters",
  "$id": "https://raw.githubusercontent.com/gazoduc/master/nextflow_schema.json",
  "description": "VCF Genotyping Pipeline",
  "type": "object",
</xsl:stylesheet>
