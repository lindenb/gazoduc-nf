<?xml version='1.0' ?>
<xsl:stylesheet
	xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
	version='1.0'
	>
<xsl:output method="text" />

<xsl:template match="/">
<xsl:apply-templates select="schema"/>
</xsl:template>

<xsl:template match="schema">
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "title":"<xsl:apply-templates select="@title"/>",
  "description":"<xsl:apply-templates select="description"/>",
  "type": "object",
  "defs": {
	<xsl:for-each select="section">
		<xsl:if test="position()&gt;1">,</xsl:if>
		<xsl:apply-templates select="."/>
	</xsl:for-each>
  },
  "allOf": [
     <xsl:for-each select="section">
	 <xsl:if test="position()&gt;1">,</xsl:if>
    	{ "$ref": "#/defs/<xsl:value-of select="@id"/>" } 
     </xsl:for-each>
  ]

}
</xsl:template>

<xsl:template match="section">
"<xsl:value-of select="@id"/>":{
	"title":"<xsl:value-of select="@title"/>",
	"type":"object",
	"required": [<xsl:for-each select="param[@required='true']">
		 <xsl:if test="position()&gt;1">,</xsl:if>
		<xsl:text>"</xsl:text>
		<xsl:value-of select="@name"/>
		<xsl:text>"</xsl:text>
		</xsl:for-each>],
	"properties": {
	<xsl:for-each select="param">
                 <xsl:if test="position()&gt;1">,</xsl:if>
		<xsl:apply-templates select="."/>
	</xsl:for-each>
	}
	
}
</xsl:template>

<xsl:template match="param">
                <xsl:text>"</xsl:text>
                <xsl:value-of select="@name"/>
                <xsl:text>":{</xsl:text>
		<xsl:for-each select="type|format|exusts|mimetype|pattern|description|fa_icon|help_text">
        	         <xsl:if test="position()&gt;1">,</xsl:if>
			<xsl:apply-templates select="."/>
		</xsl:for-each>
<xsl:text>}</xsl:text>
</xsl:template>

<xsl:template match="type|format|exusts|mimetype|pattern|description|fa_icon|help_text">
<xsl:text>
"</xsl:text>
<xsl:value-of select="name()"/>
<xsl:text>":"</xsl:text>
<xsl:value-of select="text()"/>
<xsl:text>"</xsl:text>
</xsl:template>

</xsl:stylesheet>
