<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet
 version="1.0"
 xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
 xmlns:h="http://www.w3.org/1999/xhtml"
 xmlns="http://www.w3.org/2000/svg"
 >



<!--===========================================================================-->
<xsl:template match="rec">
    <xsl:variable name="url">
        <xsl:call-template name="ucsc_anchor">
            <xsl:with-param name="build">
                    <xsl:value-of select="$build"/>
            </xsl:with-param>
            <xsl:with-param name="contig">
                    <xsl:value-of select="@contig"/>
            </xsl:with-param>
            <xsl:with-param name="start">
                    <xsl:value-of select="@start"/>
            </xsl:with-param>
            <xsl:with-param name="end">
                    <xsl:value-of select="@end"/>
            </xsl:with-param>
    </xsl:call-template>
    </xsl:variable>

    <xsl:choose>
        <xsl:when test="string-length($url) &gt; 0 and $url != '#'">
        <a>
            <xsl:attribute name="href">
                <xsl:value-of select="$url"/>
            </xsl:attribute>
            <xsl:apply-templates select="." mode="after_url"/>
        </a>
        </xsl:when>
        <xsl:otherwise>
            <xsl:apply-templates select="." mode="after_url"/>
        </xsl:otherwise>	
    </xsl:choose>

</xsl:template>
<!--===========================================================================-->

</xsl:stylesheet>
