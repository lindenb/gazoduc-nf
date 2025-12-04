<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet
 version="1.0"
 xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
 xmlns="http://www.w3.org/2000/svg"
 >
<!--

this stylesheet converts a ucsc cytoband +jvarkit bed2xml to a SVG fragment

-->


<xsl:output method="xml" encoding="UTF-8" />
<xsl:include href="circular.mod.xsl"/>
<xsl:include href="rec.mod.xsl"/>
<xsl:param name="radius_R1">1000</xsl:param>
<xsl:param name="radius_R2">990</xsl:param>

<xsl:variable name="build" select="/bed/header/dictionary/@build"/>

<xsl:template match="/">
	<g id="cytoband">
	
	<xsl:comment> BEGIN CYTOBANDS </xsl:comment>
	<xsl:apply-templates select="bed/body/rec"/>
	<xsl:comment> END CYTOBANDS </xsl:comment>
	</g>
</xsl:template>


<!--===========================================================================-->

<xsl:template match="rec" mode="after_url">

<xsl:variable name="r1" select="number($radius_R1)"/>
<xsl:variable name="r2" select="number($radius_R2)"/>


<path class="cytoband">
     <xsl:attribute name="style">
	 	<xsl:text>fill:</xsl:text>
    	 <xsl:call-template name="cytoband_color">
                <xsl:with-param name="stain">
                        <xsl:value-of select="color/text()"/>
                </xsl:with-param>
        </xsl:call-template>
		<xsl:text>;</xsl:text>
    </xsl:attribute>
    
    
    <xsl:attribute name="d">
        <xsl:call-template name="arc">
                <xsl:with-param name="f1">
                        <xsl:value-of select="number(@f1)"/>
                </xsl:with-param>
                <xsl:with-param name="f2">
                        <xsl:value-of select="number(@f2)"/>
                </xsl:with-param>
                <xsl:with-param name="outer_radius">
                        <xsl:value-of select="$r1"/>
                </xsl:with-param>
                <xsl:with-param name="inner_radius">
                        <xsl:value-of select="$r2"/>
                </xsl:with-param>
        </xsl:call-template>
	</xsl:attribute>
	<title>
		<xsl:value-of select="@contig"/>
		<xsl:text>.</xsl:text>
		<xsl:value-of select="name/text()"/>
	</title>
</path>

</xsl:template>

</xsl:stylesheet>
