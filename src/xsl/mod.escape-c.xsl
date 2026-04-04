<?xml version='1.0' encoding="UTF-8"?>
<xsl:stylesheet
        xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
        xmlns:x="http://www.ibm.com/xmlns/prod/2009/jsonx"
        version='1.0'
        >

<xsl:template name="escape-C">
<xsl:param name="s"/>
<xsl:call-template name="C-string-replace">
  <xsl:with-param name="from">"</xsl:with-param>
  <xsl:with-param name="to">\"</xsl:with-param>
  <xsl:with-param name="string">
                <xsl:call-template name="C-string-replace">
                  <xsl:with-param name="from"><xsl:text>&#10;</xsl:text></xsl:with-param>
                  <xsl:with-param name="to">\n</xsl:with-param>
                  <xsl:with-param name="string">
                        <xsl:call-template name="C-string-replace">
                          <xsl:with-param name="from">\</xsl:with-param>
                          <xsl:with-param name="to">\\</xsl:with-param>
                          <xsl:with-param name="string" select="$s"/>
                        </xsl:call-template>
                  </xsl:with-param>
                </xsl:call-template>
  </xsl:with-param>
</xsl:call-template>
</xsl:template>


<xsl:template name="C-string-replace">
    <xsl:param name="string"/>
    <xsl:param name="from"/>
    <xsl:param name="to"/>
    
    <xsl:if test="string-length($from)=0"><xsl:message terminate="yes">BOUM:'<xsl:value-of select="$from"/>' vs '<xsl:value-of select="$to"/>'</xsl:message></xsl:if>
    
    <xsl:choose>
      <xsl:when test="contains($string,$from)">
        <xsl:value-of select="substring-before($string,$from)"/>
        <xsl:value-of select="$to"/>
        <xsl:call-template name="json-string-replace">
          <xsl:with-param name="string" select="substring-after($string,$from)"/>
          <xsl:with-param name="from" select="$from"/>
          <xsl:with-param name="to" select="$to"/>
        </xsl:call-template>
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="$string"/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

</xsl:stylesheet>
