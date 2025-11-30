<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet 
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:xi="http://www.w3.org/2001/XInclude"
    version="1.0"
    >
  <xsl:output method="html" encoding="UTF-8"/>
  
  <xsl:param name="title"></xsl:param>

  <xsl:template match="/">
    <xsl:apply-templates/>
  </xsl:template>

  <xsl:template match="title[name(..)='head']">
    <title><xsl:value-of select="$title"/></title>
  </xsl:template>
   
  <xsl:template match="body[name(..)='html']">
  <body>
  <h1><xsl:value-of select="$title"/></h1>
  <div>
   <xi:include href="header.xml"/>
  </div>
  <xsl:apply-templates/>
  <div id="vcf2table">
   <xi:include href="footer.xml"/>
  </div>
  </body>
  </xsl:template>
  
  <xsl:template match="style[name(..)='head']">
        <xsl:copy>
          <xsl:copy-of select="@*"/>
          <xsl:apply-templates/>
	  <xsl:text>
table.minimalistBlack { border: 1px solid #1F1F1F; text-align: left; border-collapse: collapse; }
table.minimalistBlack td, table.minimalistBlack th { border: 1px solid #1F1F1F; padding: 5px 2px; }
table.minimalistBlack tbody td { font-size: 13px; }
table.minimalistBlack thead { background: #CFCFCF; background: -moz-linear-gradient(top, #dbdbdb 0%, #d3d3d3 66%, #CFCFCF 100%);background: -webkit-linear-gradient(top, #dbdbdb 0%, #d3d3d3 66%, #CFCFCF 100%); background: linear-gradient(to bottom, #dbdbdb 0%, #d3d3d3 66%, #CFCFCF 100%); border-bottom: 2px solid #000000; }
table.minimalistBlack thead th { font-size: 15px; font-weight: bold; color: #000000; text-align: left; }
table.minimalistBlack tfoot td { font-size: 14px; } 
</xsl:text>
        </xsl:copy>
  </xsl:template>


  <xsl:template match="@*|node()">
        <xsl:copy>
          <xsl:copy-of select="@*"/>
          <xsl:apply-templates/>
        </xsl:copy>
  </xsl:template>

</xsl:stylesheet>