<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet
 version="1.0"
 xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
 xmlns="http://www.w3.org/2000/svg"
 exclude-result-prefixes=""
 >
<!--

this stylesheet converts a ucsc cytoband +jvarkit bed2xml to a SVG fragment

-->


<xsl:output method="xml" encoding="UTF-8"  omit-xml-declaration="yes" indent="yes" />
<xsl:include href="circular.mod.xsl"/>
<xsl:include href="rec.mod.xsl"/>
<xsl:param name="radius_R1">1000</xsl:param>
<xsl:param name="radius_R2">990</xsl:param>
<xsl:param name="title"></xsl:param>
<xsl:param name="style"></xsl:param>
<xsl:param name="class"></xsl:param>

<xsl:variable name="build" select="/bed/header/dictionary/@build"/>

<!--===========================================================================-->

<xsl:template match="/">
	<xsl:variable name="r1" select="number($radius_R1)"/>
    <xsl:variable name="r2" select="number($radius_R2)"/>
    <xsl:variable name="height" select="($r1 - $r2)"/>
    <xsl:variable name="mid" select="$r1 - ($height div 2)"/>

	<g>
        <xsl:comment> BEGIN INDEXCOV </xsl:comment>


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

        <xsl:apply-templates select="bed/body/rec"/>

        


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

        <circle cx="0" cy="0" class="MID">
            <xsl:attribute name="r">
                 <xsl:value-of select="$mid"/>
            </xsl:attribute>
        </circle>

        <circle cx="0" cy="0" class="DEL">
            <xsl:attribute name="r">
                 <xsl:value-of select="$mid - ($height div 4.0) "/>
            </xsl:attribute>
        </circle>

         <circle cx="0" cy="0" class="DUP">
            <xsl:attribute name="r">
                 <xsl:value-of select="$mid + ($height div 4.0) "/>
            </xsl:attribute>
        </circle>


        <xsl:comment> END INDEXCOV </xsl:comment>
	</g>
</xsl:template>



<!--===========================================================================-->


<xsl:template match="rec" mode="after_url">
<xsl:variable name="norm0" select=" number(norm/text())"/>
<xsl:variable name="norm">
    <xsl:choose>
        <xsl:when test="$norm0 &lt; 0.0"><xsl:value-of select="number(0.0)"/></xsl:when>
        <xsl:when test="$norm0 &gt; 2.0"><xsl:value-of select="number(2.0)"/></xsl:when>
        <xsl:otherwise><xsl:value-of select="$norm0"/></xsl:otherwise>
    </xsl:choose>
</xsl:variable>
<xsl:variable name="rA" select="number($radius_R1)"/>
<xsl:variable name="rB" select="number($radius_R2)"/>
<xsl:variable name="height" select="($rA - $rB)"/>
<xsl:variable name="mid" select="$rA - ($height div 2)"/>


<xsl:variable name="r1">
    <xsl:choose>
        <xsl:when test="$norm &gt; 1.0">
            <xsl:value-of select="$mid + ($norm - 1.0) * ( $height  div 2.0) "/>
        </xsl:when>
        <xsl:otherwise><xsl:value-of select="$mid"/></xsl:otherwise>
    </xsl:choose>
</xsl:variable>

<xsl:variable name="r2">
    <xsl:choose>
        <xsl:when test="$norm &lt; 1.0">
            <xsl:value-of select="$mid - norm * ( $height  div 2.0) "/>
        </xsl:when>
        <xsl:otherwise><xsl:value-of select="$mid"/></xsl:otherwise>
    </xsl:choose>
</xsl:variable>


<path>

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
            <xsl:value-of select="norm"/>
        </title>
</path>
</xsl:template>

</xsl:stylesheet>
