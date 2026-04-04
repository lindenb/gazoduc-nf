<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet
        version="1.0"
        xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
 >
<xsl:include href="mod.escape-c.xsl"/>
<xsl:output method="text" encoding="UTF-8" />

<!-- ========================================================================================== -->
<!-- ========================================================================================== -->
<!-- ========================================================================================== -->

<xsl:template match="/">
<xsl:apply-templates select="schema"/>
</xsl:template>

<!-- ========================================================================================== -->
<!-- ========================================================================================== -->
<!-- ========================================================================================== -->

<xsl:template match="schema">
{
   "$id" :   <xsl:choose>
	<xsl:when test="id">
		<xsl:apply-templates select="id" mode="quote"/>
	</xsl:when>
	<xsl:otherwise>
		<xsl:text>"https://gitlab.univ-nantes.fr/pierre.lindenbaum/gazoduc-nf/-/tree/master/workflows/</xsl:text>
    <xsl:value-of select="generate-id(.)"/>
    <xsl:text>"</xsl:text>
	</xsl:otherwise>
  </xsl:choose>,
  "title" :   <xsl:choose>
	<xsl:when test="@name">
		<xsl:apply-templates select="@name" mode="quote"/>
	</xsl:when>
	<xsl:otherwise>
		<xsl:text>"undefined workflow"</xsl:text>
	</xsl:otherwise>
  </xsl:choose>,
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "description": <xsl:choose>
	<xsl:when test="description">
		<xsl:apply-templates select="description" mode="quote"/>
	</xsl:when>
	<xsl:otherwise>
		<xsl:text>"no description"</xsl:text>
	</xsl:otherwise>
  </xsl:choose>,
  "type": "object",
  "$defs": {
    <xsl:for-each select="params/section">
    <xsl:if test="position() &gt; 1">
      <xsl:text>,</xsl:text>
    </xsl:if>

    <xsl:text>"</xsl:text>
    <xsl:value-of select="generate-id(.)"/>
    <xsl:text>" : {</xsl:text>
    <xsl:apply-templates select="."/>
    <xsl:text> }</xsl:text>
    </xsl:for-each>
    },
"allOf": [

     <xsl:for-each select="params/section">
      <xsl:if test="position() &gt; 1">
      <xsl:text>,</xsl:text>
    </xsl:if>
    <xsl:text>{
      "$ref": "#/$defs/</xsl:text>
    <xsl:value-of select="generate-id(.)"/>
    <xsl:text>"}</xsl:text>
     </xsl:for-each>
  ]
}
</xsl:template>


<!-- ========================================================================================== -->
<!-- ========================================================================================== -->
<!-- ========================================================================================== -->


<xsl:template match="section">
 "title" :  <xsl:choose>
	<xsl:when test="@name">
		<xsl:apply-templates select="@name" mode="quote"/>
	</xsl:when>
	<xsl:otherwise>
		<xsl:text>"undefined section"</xsl:text>
	</xsl:otherwise>
  </xsl:choose>,
"description": <xsl:choose>
	<xsl:when test="description">
		<xsl:apply-templates select="description" mode="quote"/>
	</xsl:when>
	<xsl:otherwise>
		<xsl:text>"no description"</xsl:text>
	</xsl:otherwise>
  </xsl:choose>,
"type": "object",
"fa_icon":  <xsl:choose>
	<xsl:when test="fa_icon">
		<xsl:apply-templates select="fa_icon" mode="quote"/>
	</xsl:when>
	<xsl:otherwise>
		<xsl:text>"fas fa-terminal"</xsl:text>
	</xsl:otherwise>
  </xsl:choose>,
<xsl:if test="help_text">
"help_text": 	<xsl:apply-templates select="help_text" mode="quote"/>,
</xsl:if>
"required":[
  <xsl:for-each select="param[@required='true']">
      <xsl:if test="position() &gt; 1">
        <xsl:text>,</xsl:text>
      </xsl:if>
      <xsl:apply-templates select="@name" mode="quote"/>
  </xsl:for-each>
],
"properties":{
   <xsl:for-each select="param">
      <xsl:if test="position() &gt; 1">
        <xsl:text>,</xsl:text>
      </xsl:if>
      <xsl:apply-templates select="."/>
  </xsl:for-each>
}
</xsl:template>

<xsl:template match="param">
<xsl:apply-templates select="@name" mode="quote"/>
<xsl:text>: {
</xsl:text>
"type":<xsl:choose>
	<xsl:when test="type">
		<xsl:apply-templates select="type" mode="quote"/>
	</xsl:when>
  
	<xsl:otherwise>
		<xsl:text>"string"</xsl:text>
	</xsl:otherwise>
  </xsl:choose>,
"description": <xsl:choose>
	<xsl:when test="description">
		<xsl:apply-templates select="description" mode="quote"/>
	</xsl:when>
	<xsl:otherwise>
		<xsl:apply-templates select="@name" mode="quote"/>
	</xsl:otherwise>
  </xsl:choose>
<xsl:choose>
  <xsl:when test="default">
    <xsl:text>,"default":</xsl:text>
    <xsl:apply-templates select="default" mode="quote"/>
  </xsl:when>
  <xsl:when test="type='boolean'">
    <xsl:text>,"default":false</xsl:text>
  </xsl:when>
</xsl:choose>
<xsl:if test="format">
,"format": 	<xsl:apply-templates select="format" mode="quote"/>
</xsl:if>
<xsl:if test="pattern">
,"pattern": 	<xsl:apply-templates select="pattern" mode="quote"/>
</xsl:if>
<xsl:if test="help_text">
,"help_text": 	<xsl:apply-templates select="help_text" mode="quote"/>
</xsl:if>

<xsl:choose>
  <xsl:when test="fa_icon">
    <xsl:text>,"fa_icon":</xsl:text>
    <xsl:apply-templates select="fa_icon" mode="quote"/>
  </xsl:when>
  <xsl:when test="type='boolean'">
    <xsl:text>,"fa_icon" : "fas fa-question-circle"</xsl:text>
  </xsl:when>
</xsl:choose>

<xsl:text>
}</xsl:text>
</xsl:template>

<!-- ========================================================================================== -->
<!-- ========================================================================================== -->
<!-- ========================================================================================== -->

<xsl:template match="*|@*" mode="quote">
<xsl:text>"</xsl:text>
<xsl:call-template name="escape-C">
  <xsl:with-param name="s" select="."/>
</xsl:call-template>
<xsl:text>"</xsl:text>
</xsl:template>

<!-- ========================================================================================== -->
<!-- ========================================================================================== -->
<!-- ========================================================================================== -->

</xsl:stylesheet>
