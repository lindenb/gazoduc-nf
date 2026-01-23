<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet
	version="1.0"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
	xmlns:foaf="http://xmlns.com/foaf/0.1/"
	xmlns:dc="http://purl.org/dc/elements/1.1/"
	xmlns:u="http://www.x"
 >



<xsl:output method="xml" encoding="UTF-8" />




<!--===========================================================================-->

<xsl:template match="/">
	<rdf:RDF>
	<xsl:apply-templates select="//gene[attributes/attribute[@key='gene_id']]"/>
	</rdf:RDF>
</xsl:template>



<!--===========================================================================-->

<xsl:template match="gene" >
<u:gene>
<xsl:attribute name="rdf:about">
	
	
	<xsl:text>u:</xsl:text>
	<xsl:value-of select="attributes/attribute[@key='gene_id']/text()"/>
</xsl:attribute>
<xsl:apply-templates select="." mode="coord"/>

<xsl:for-each select="//transcript">
	<u:has_transcript>
		<xsl:apply-templates select="."/>
	</u:has_transcript>
</xsl:for-each>

</u:gene>
</xsl:template>


<!--===========================================================================-->

<xsl:template match="transcript" >
<u:gene>
<xsl:attribute name="rdf:about">
	
	
	<xsl:text>u:</xsl:text>
	<xsl:value-of select="attributes/attribute[@key='transcript_id']/text()"/>
</xsl:attribute>
<xsl:apply-templates select="." mode="coord"/>

<xsl:for-each select="//exon">
	<u:has_exon>
	</u:has_exon>
</xsl:for-each>

</u:gene>
</xsl:template>

<xsl:template match="*" mode="coord" >
<u:chrom><xsl:value-of select="@chrom"/></u:chrom>
<u:start xsdType=""><xsl:value-of select="@start"/></u:start>
<u:end xsdType=""><xsl:value-of select="@end"/></u:end>

</xsl:template>


</xsl:stylesheet>
