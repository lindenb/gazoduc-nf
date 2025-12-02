<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet
 version="1.0"
 xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
 xmlns="http://www.w3.org/2000/svg"
 >
<!--

this stylesheet converts a faix+jvarkit bed2xml to a SVG fragment

-->


<xsl:output method="xml" encoding="UTF-8" />
<xsl:include href="circular.mod.xsl"/>
<xsl:param name="radius">1000</xsl:param>

<xsl:template match="/">
<xsl:apply-template select="bed/body"/>
</xsl:template>

<xsl:template match="body">
<g>
<xsl:apply-template select="rec"/>
</g>
</xsl:template>


<xsl:template match="header">
<g>
<xsl:apply-templates select="dictionary/sequence"/>
</g>
</xsl:template>

</xsl:stylesheet>
