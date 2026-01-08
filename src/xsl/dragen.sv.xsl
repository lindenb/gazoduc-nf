<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet
 version="1.0"
 xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
 xmlns="http://www.w3.org/2000/svg"
 >
<!--

this stylesheet converts a bed2xml to a SVG fragment

-->


<xsl:output method="xml" encoding="UTF-8" omit-xml-declaration="yes" />
<xsl:include href="circular.mod.xsl"/>
<xsl:include href="rec.mod.xsl"/>
<xsl:param name="radius_R1">1000</xsl:param>
<xsl:param name="radius_R2">990</xsl:param>
<xsl:param name="style"></xsl:param>
<xsl:param name="title"></xsl:param>

<xsl:variable name="build" select="/bed/header/dictionary/@build"/>

<!--===========================================================================-->

<xsl:template match="/">
	<g>
		<xsl:attribute name="id">
		<xsl:text>CNV</xsl:text>
		<xsl:value-of select="$title"/>
		</xsl:attribute>
	
	<xsl:call-template name="donut">
		<xsl:with-param name="outer_radius">
			<xsl:value-of select="$radius_R1"/>
		</xsl:with-param>
		<xsl:with-param name="inner_radius">
			<xsl:value-of select="$radius_R2"/>
		</xsl:with-param>
		<xsl:with-param name="style">
			<xsl:value-of select="$style"/>
		</xsl:with-param>
		<xsl:with-param name="title">
			<xsl:value-of select="$title"/>
		</xsl:with-param>
	</xsl:call-template>


	<circle cx="0" cy="0" r="{$radius_R1}" style="fill:none;stroke:lightgray;stroke-dasharray:0 4 0;"/>
	<xsl:comment> BEGIN CNV </xsl:comment>
	<xsl:apply-templates select="bed/body/contig/row/rec"/>
	<xsl:comment> END CNV </xsl:comment>
	<text x="0">
		<xsl:attribute name="y">
		<xsl:value-of select="($radius_R1 + $radius_R2) div 2.0"/>			
		</xsl:attribute>
		<xsl:attribute name="style">
			<xsl:text>dominant-baseline:middle;text-anchor:middle;font-size:</xsl:text>
			<xsl:value-of select="($radius_R1 - $radius_R2) * 0.8"/>		
			<xsl:text>px;</xsl:text>
		</xsl:attribute>
		<xsl:value-of select="$title"/>
	</text>
	</g>
</xsl:template>




<!--===========================================================================-->


<xsl:template match="rec" mode="after_url">
<xsl:variable name="nrows" select="number(../../@rows)"/>
<xsl:variable name="dr" select="(number($radius_R1) - number($radius_R2)) div $nrows"/>
<xsl:variable name="r1" select="number($radius_R1) - ( number(@y) * $dr)"/>
<xsl:variable name="r2" select="$r1 - ($dr * 0.95)"/>

<path >
        <xsl:attribute name="class">
                <xsl:text>cnv </xsl:text>
		<xsl:value-of select="svtype"/>
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
		
       		<title><xsl:value-of select="svtype"/></title>
</path>
</xsl:template>

</xsl:stylesheet>
