<?xml version='1.0' ?>
<xsl:stylesheet
	xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
	version='1.0'
	>
<xsl:output method="xml" ident="yes"/>

<xsl:template match="/">
<xsl:apply-templates select="schema"/>
</xsl:template>


<xsl:template match="macro[@id='help']">
	<param name="help">
                        <type>boolean</type>
                        <description>display help</description>
                        <help_text>display help</help_text>
	</param>
</xsl:template>

<xsl:template match="macro[@id='outdir']">
 <param name="outdir" required="true">
    <type>string</type>
    <format>directory-path</format>
    <exists>true</exists>
    <mimetype>text/plain</mimetype>
    <description>The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure</description>
    <fa_icon>far fa-file-code</fa_icon>
 </param>
</xsl:template>

<xsl:template match="macro[@id='fasta']">
 <param name="fasta">
    <xsl:if test="@required='true'"><xsl:attribute name="required">true</xsl:attribute></xsl:if>
    <type>string</type>
    <format>file-path</format>
    <exists>true</exists>
    <mimetype>text/plain</mimetype>
    <pattern>^\\S+\\.fn?a(sta)?(\\.gz)?$</pattern>
    <description>Path to FASTA genome file</description>
    <fa_icon>far fa-file-code</fa_icon>
 </param>
</xsl:template>

<xsl:template match="macro[@id='prefix']">
  <param name="prefix" required="true">
     <type>string</type>
     <mimetype>text/plain</mimetype>
     <description>file prefix for output files</description>
  </param>
</xsl:template>

<xsl:template match="macro[@id='mosdepth_version']">
  <param>
     <xsl:attribute name="name"><xsl:value-of select="@id"/></xsl:attribute>
     <type>string</type>
     <mimetype>text/plain</mimetype>
     <description>mosdepth version</description>
  </param>
</xsl:template>

<xsl:template match="macro[@id='somalier_version']">
  <param>
     <xsl:attribute name="name"><xsl:value-of select="@id"/></xsl:attribute>
     <type>string</type>
     <mimetype>text/plain</mimetype>
     <description>somalier version</description>
  </param>
</xsl:template>


<xsl:template match="macro[@id='mapq']">
  <param>
     <xsl:attribute name="name"><xsl:value-of select="@id"/></xsl:attribute>
     <type>integer</type>
     <description>min mapping quality MAPQ</description>
  </param>
</xsl:template>



<xsl:template match="*">
<xsl:copy>
<xsl:copy-of select="@*"/>
<xsl:apply-templates select="*|text()"/>
</xsl:copy>
</xsl:template>

<xsl:template match="text()">
<xsl:copy>
<xsl:apply-templates select="*|text()"/>
</xsl:copy>
</xsl:template>

</xsl:stylesheet>
