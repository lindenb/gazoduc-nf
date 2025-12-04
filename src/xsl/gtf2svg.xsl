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

<!--===========================================================================-->

<xsl:template match="/">
	<g id="cytoband">
	<xsl:comment> BEGIN GTF </xsl:comment>
	<xsl:apply-templates select="bed/body/contig/row/rec"/>
	<xsl:comment> END GTF </xsl:comment>
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
                <xsl:text>gene</xsl:text>
						<xsl:text> </xsl:text>
				<xsl:value-of select="gene_biotype"/>
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
		
		<xsl:if test='gene_name/text()!="." and gene_name/text()!=""'>
       		<title><xsl:value-of select="gene_name"/></title>
		</xsl:if>
</path>
</xsl:template>

</xsl:stylesheet>
