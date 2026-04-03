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

<xsl:template match="gene[attributes/attribute[@key='gene_id']]" >
<u:gene>
<xsl:attribute name="rdf:about">
	<xsl:text>u:</xsl:text>
	<xsl:value-of select="attributes/attribute[@key='gene_id']/text()"/>
</xsl:attribute>

<xsl:for-each select="//transcript">

</xsl:for-each>

</u:gene>
</xsl:template>

</xsl:stylesheet>
