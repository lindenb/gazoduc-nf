<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet
	version="1.0"
	xmlns:obo="http://purl.obolibrary.org/obo/"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
	xmlns:foaf="http://xmlns.com/foaf/0.1/"
	xmlns:dc="http://purl.org/dc/elements/1.1/"
	xmlns:u="https://umr1087.univ-nantes.fr/rdf/"
 >



<xsl:output method="xml" encoding="UTF-8" indent="yes"/>

<xsl:include href="mod.rec.biotype.xslt.xsl"/><!-- generated with SPARQL+XSLT-->


<!--===========================================================================-->

<xsl:template match="/">
	<rdf:RDF>
	<xsl:apply-templates select="gtf"/>
	</rdf:RDF>
</xsl:template>

<!--===========================================================================-->

<xsl:template match="gtf">
	<xsl:apply-templates select="dictionary"/>
	<xsl:apply-templates select="//gene[attributes/attribute[@key='gene_id']]"/>
</xsl:template>



<!--===========================================================================-->
<xsl:template match="dictionary">
	<obo:SO_0001505>
		<xsl:attribute name="rdf:about">
			<xsl:text>https://umr1087.univ-nantes.fr/rdf/build/</xsl:text>
			<xsl:value-of select="../@build"/>
		</xsl:attribute>
		<dc:title><xsl:value-of select="@build"/></dc:title>
	</obo:SO_0001505>
	
	<xsl:apply-templates select="sequence"/>
</xsl:template>

<!--===========================================================================-->

<xsl:template match="sequence" >
<obo:SO_0000340>
<xsl:attribute name="rdf:about">
	<xsl:text>https://umr1087.univ-nantes.fr/rdf/build/</xsl:text>
	<xsl:value-of select="../@build"/>
	<xsl:text>/</xsl:text>
	<xsl:value-of select="@name"/>
</xsl:attribute>
<dc:title><xsl:value-of select="@name"/></dc:title>
<u:length  rdf:datatype="http://www.w3.org/2001/XMLSchema#int"><xsl:value-of select="@name"/></u:length>
</obo:SO_0000340>
</xsl:template>


<!--===========================================================================-->

<xsl:template match="gene" >
<xsl:variable name="type">
	<xsl:value-of select="attributes/attribute[@key='gene_biotype']"/>
</xsl:variable>

<xsl:variable name="element_name">
<xsl:call-template name="biotype2qname">
  <xsl:with-param name="type" select="$type"/>
  <xsl:with-param name="default" ><xsl:text>obo:SO_0000704</xsl:text></xsl:with-param>
</xsl:call-template>
</xsl:variable>


<xsl:element name="{$element_name}">
		<xsl:apply-templates select="." mode="gene_content"/>
</xsl:element>

</xsl:template>


<xsl:template match="gene" mode="gene_content">

<xsl:attribute name="rdf:about">
	
	
	<xsl:text>http://rdf.ebi.ac.uk/resource/ensembl/</xsl:text>
	<xsl:value-of select="attributes/attribute[@key='gene_id']/text()"/>
</xsl:attribute>
<xsl:apply-templates select="." mode="coord"/>
<xsl:apply-templates select="attributes"/>

<xsl:for-each select=".//transcript">
	<u:has_transcript>
		<xsl:apply-templates select="."/>
	</u:has_transcript>
</xsl:for-each>


</xsl:template>


<!--===========================================================================-->

<xsl:template match="transcript" >
<obo:SO_0000673>
<xsl:attribute name="rdf:about">
	<xsl:text>http://rdf.ebi.ac.uk/resource/ensembl/</xsl:text>
	<xsl:value-of select="attributes/attribute[@key='transcript_id']/text()"/>
</xsl:attribute>
<xsl:apply-templates select="." mode="coord"/>
<xsl:apply-templates select="attributes"/>

<xsl:for-each select=".//exon">
	<u:has_exon>
		<xsl:apply-templates select="."/>
	</u:has_exon>
</xsl:for-each>

<xsl:for-each select=".//cds">
	<u:has_cds>
		<xsl:apply-templates select="."/>
	</u:has_cds>
</xsl:for-each>


</obo:SO_0000673>
</xsl:template>

<!--===========================================================================-->

<xsl:template match="exon" >
<obo:SO_0000147>
<xsl:apply-templates select="." mode="coord"/>
</obo:SO_0000147>
</xsl:template>

<!--===========================================================================-->


<xsl:template match="cds" >
<obo:SO_0000316>
<xsl:apply-templates select="." mode="coord"/>
</obo:SO_0000316>
</xsl:template>

<!--===========================================================================-->

<xsl:template match="attributes" >
<xsl:apply-templates select="attribute"/>
</xsl:template>


<xsl:template match="attribute" >
	<xsl:choose>
		<xsl:when test="@key='gene_id'"></xsl:when>
		<xsl:when test="@key='transcript_id'"></xsl:when>
		<xsl:when test="@key='exon_id'"></xsl:when>
		<xsl:when test="@key='exon_version'"></xsl:when>
		<xsl:when test="@key='gene_version'"></xsl:when>
		<xsl:when test="@key='gene_source'"></xsl:when>
		<xsl:when test="@key='transcript_source'"></xsl:when>
		<xsl:when test="@key='tag'"></xsl:when>
		<xsl:when test="@key='transcript_support_level'"></xsl:when>
		<xsl:otherwise>
			<xsl:element name="u:{@key}">
				<xsl:value-of select="./text()"/>
			</xsl:element>
		</xsl:otherwise>
	</xsl:choose>
</xsl:template>

<!--===========================================================================-->



<xsl:template match="*" mode="coord" >
<u:has_interval>
<u:Interval>
<u:chrom>
<xsl:attribute name="rdf:resource">
	<xsl:text>https://umr1087.univ-nantes.fr/rdf/build/</xsl:text>
	<xsl:value-of select="/gtf/dictionary/@build"/>
	<xsl:text>/</xsl:text>
	<xsl:value-of select="@chrom"/>
</xsl:attribute>
</u:chrom>
<u:start rdf:datatype="http://www.w3.org/2001/XMLSchema#int"><xsl:value-of select="@start"/></u:start>
<u:end rdf:datatype="http://www.w3.org/2001/XMLSchema#int"><xsl:value-of select="@end"/></u:end>
</u:Interval>
</u:has_interval>
</xsl:template>


</xsl:stylesheet>
