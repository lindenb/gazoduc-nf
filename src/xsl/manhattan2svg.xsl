<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet
 version="1.0"
 xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
 xmlns="http://www.w3.org/2000/svg"
 >
<!--

this stylesheet converts a ucsc cytoband +jvarkit bed2xml to a SVG fragment

-->


<xsl:output method="xml" encoding="UTF-8"  omit-xml-declaration="yes" indent="yes" />
<xsl:include href="circular.mod.xsl"/>
<xsl:include href="rec.mod.xsl"/>
<xsl:param name="radius_R1">900</xsl:param>
<xsl:param name="radius_R2">100</xsl:param>
<xsl:param name="title"></xsl:param>
<xsl:param name="style"></xsl:param>
<xsl:param name="class"></xsl:param>
<xsl:param name="min_value"></xsl:param>
<xsl:param name="max_value"></xsl:param>

<xsl:variable name="build" select="/bed/header/dictionary/@build"/>

<!--===========================================================================-->
<xsl:variable name="MAX_VALUE">
    <xsl:choose>
        <xsl:when test="string-length($max_value) &gt; 0">
            <xsl:value-of select="number($max_value)"/>
        </xsl:when>
        <xsl:when test="/bed/body/rec/value">
            <xsl:for-each select="/bed/body/rec/value">
                <xsl:sort select="." data-type="number" order="descending"/>
                <xsl:if test="position() = 1">
                    <xsl:value-of select="number(.)"/>
                </xsl:if>
            </xsl:for-each>
        </xsl:when>
        <xsl:otherwise>
            <xsl:value-of select="number(1.0)"/>
        </xsl:otherwise>
    </xsl:choose>
</xsl:variable>

<!--===========================================================================-->
<xsl:variable name="MIN_VALUE">
    <xsl:choose>
        <xsl:when test="string-length($min_value) &gt; 0">
            <xsl:value-of select="number($min_value)"/>
        </xsl:when>
        <xsl:when test="/bed/body/rec/value">
            <xsl:for-each select="/bed/body/rec/value">
                <xsl:sort select="." data-type="number" order="ascending"/>
                <xsl:if test="position() = 1">
                    <xsl:value-of select="number(.)"/>
                </xsl:if>
            </xsl:for-each>
        </xsl:when>
        <xsl:otherwise>
            <xsl:value-of select="number(0)"/>
        </xsl:otherwise>
    </xsl:choose>
</xsl:variable>

<!--===========================================================================-->

<xsl:template match="/">
	

	<g>
        <xsl:comment> BEGIN MANHATTAN </xsl:comment>


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
            <xsl:with-param name="class">
                <xsl:value-of select="$class"/>
            </xsl:with-param>
            <xsl:with-param name="title">
                <xsl:value-of select="$title"/>
            </xsl:with-param>
        </xsl:call-template>

        <xsl:apply-templates select="bed/body/rec[value]"/>



        <xsl:comment> END MANHATTAN </xsl:comment>
	</g>
</xsl:template>



<!--===========================================================================-->


<xsl:template match="rec" mode="after_url">
<xsl:variable name="norm" select=" (number(value/text())  - number($MIN_VALUE)) div (number($MAX_VALUE) - number($MIN_VALUE))"/>
<xsl:variable name="r" select="number($radius_R2) + $norm * ( number($radius_R1) -  number($radius_R2) )"/>



<circle>

        <xsl:if test="string-length(normalize-space(class)) &gt; 0">
            <xsl:attribute name="class">
                <xsl:value-of select="normalize-space(class)"/>
            </xsl:attribute>
        </xsl:if>
         <xsl:if test="string-length(normalize-space(style)) &gt; 0">
            <xsl:attribute name="style">
                <xsl:value-of select="normalize-space(style)"/>
            </xsl:attribute>
        </xsl:if>
        

    <xsl:attribute name="cx">
                <xsl:call-template name="toX">
                    <xsl:with-param name="f">
                        <xsl:value-of select="(number(@f1) + number(@f2)) div 2.0"/>
                    </xsl:with-param>
                    <xsl:with-param name="r">
                        <xsl:value-of select="$r"/>
                    </xsl:with-param>
                </xsl:call-template>
    </xsl:attribute>

    <xsl:attribute name="cy">
                <xsl:call-template name="toY">
                    <xsl:with-param name="f">
                        <xsl:value-of select="(number(@f1) + number(@f2)) div 2.0"/>
                    </xsl:with-param>
                    <xsl:with-param name="r">
                        <xsl:value-of select="$r"/>
                    </xsl:with-param>
                </xsl:call-template>
    </xsl:attribute>

    <xsl:attribute name="r">
        <xsl:choose>
            <xsl:when test="radius/text() and number(radius/text()) &gt;= 0">
                <xsl:value-of select="radius/text()"/>
            </xsl:when>
            <xsl:otherwise>
                <xsl:value-of select="number(5)"/>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:attribute>

    <xsl:if test="string-length(normalize-space(title)) &gt; 0">
        <title>
            <xsl:value-of select="title/text()"/>
        </title>
	</xsl:if>
</circle>
</xsl:template>

</xsl:stylesheet>
