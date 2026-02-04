<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE stylesheet [
<!ENTITY u  "https://umr1087.univ-nantes.fr/rdf/" >
<!ENTITY obo  "http://purl.obolibrary.org/obo/" >
<!ENTITY xsd  "http://www.w3.org/2001/XMLSchema#" >
]>
<xsl:stylesheet
	version="1.0"
	xmlns:obo="&obo;"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
	xmlns:foaf="http://xmlns.com/foaf/0.1/"
	xmlns:owl ="http://www.w3.org/2002/07/owl#"
	xmlns:dc="http://purl.org/dc/elements/1.1/"
	xmlns:u="&u;"
 >



<xsl:output method="xml" encoding="UTF-8" indent="yes"/>

<xsl:include href="mod.rec.biotype.xslt.xsl"/><!-- generated with SPARQL+XSLT-->


<!--===========================================================================-->

<xsl:template match="/">
	<rdf:RDF>
	<xsl:apply-templates select="gtf"/>

	<owl:Class rdf:about="&u;Build">
		<rdfs:label>build</rdfs:label>
		<owl:equivalentClass rdf:resource="&obo;SO_0001505"/>
	</owl:Class>

	<owl:Class rdf:about="&u;Chromosome">
		<rdfs:label>chromosome</rdfs:label>
		<owl:equivalentClass rdf:resource="&obo;SO_0000340"/>
	</owl:Class>

	<owl:Class rdf:about="&u;Gene">
		<rdfs:label>gene</rdfs:label>
		<owl:equivalentClass rdf:resource="&obo;SO_0000704"/>
	</owl:Class>

	<owl:Class rdf:about="&u;Transcript">
		<rdfs:label>transcript</rdfs:label>
		<owl:equivalentClass rdf:resource="&obo;SO_0000673"/>
	</owl:Class>

	<owl:Class rdf:about="&u;Exon">
		<rdfs:label>exon</rdfs:label>
		<owl:equivalentClass rdf:resource="&obo;S0_0000147"/>
	</owl:Class>

	<owl:Class rdf:about="&u;CDS">
		<rdfs:label>CDS</rdfs:label>
		<owl:equivalentClass rdf:resource="&obo;SO_0000316"/>
	</owl:Class>


	</rdf:RDF>
</xsl:template>

<!--===========================================================================-->

<xsl:template match="gtf">
	<xsl:apply-templates select="dictionary"/>
	<xsl:apply-templates select="gene|transcript|exon|cds"/>
</xsl:template>



<!--===========================================================================-->
<xsl:template match="dictionary">
	<u:Build>
		<xsl:attribute name="rdf:about">
			<xsl:text>&u;build</xsl:text>
			<xsl:value-of select="../@build"/>
		</xsl:attribute>
		<dc:title><xsl:value-of select="@build"/></dc:title>
	</u:Build>
	
	<xsl:apply-templates select="sequence"/>
</xsl:template>

<!--===========================================================================-->

<xsl:template match="sequence" >
<u:Chromosome>
<xsl:attribute name="rdf:about">
	<xsl:text>&u;build/</xsl:text>
	<xsl:value-of select="../@build"/>
	<xsl:text>/</xsl:text>
	<xsl:value-of select="@name"/>
</xsl:attribute>
<dc:title><xsl:value-of select="@name"/></dc:title>
<u:length  rdf:datatype="&xsd;int"><xsl:value-of select="@length"/></u:length>
</u:Chromosome>
</xsl:template>


<!--===========================================================================-->

<xsl:template match="gene|cds|transcript|exon" >

<xs:variable name="root_name">
	<xsl:choose>
		<xsl:when test="name() = 'gene' ">u:Gene</xsl:when>
		<xsl:when test="name() = 'cds' ">u:CDS</xsl:when>
		<xsl:when test="name() = 'transcript' ">u:Transcript</xsl:when>
		<xsl:when test="name() = 'exon' ">u:Exon</xsl:when>
		<xsl:otherwise>
			<xsl:message terminate="yes">
				<xsl:text>unknown GTF type</xsl:text>
				<xl:value-of select="name(.)"/>
			</xsl:message>
		</xsl:otherwise>
	</xsl:choose>
</xs:variable>


<xsl:element name="{$root_name}">

	<!-- add rdf:about for gene and transcript -->
	<xsl:choose>
		<xsl:when test="name() = 'gene' and attributes/attribute[@key='gene_id'] ">
			<xsl:attribute name="rdf:about">
				<xsl:text>http://rdf.ebi.ac.uk/resource/ensembl/</xsl:text>
				<xsl:value-of select="substring-before(attributes/attribute[@key='gene_id']/text(),'.')"/>
			</xsl:attribute>
		</xsl:when>
		<xsl:when test="name() = 'transcript' and attributes/attribute[@key='transcript_id'] ">
			<xsl:attribute name="rdf:about">
				<xsl:text>http://rdf.ebi.ac.uk/resource/ensembl/</xsl:text>
				<xsl:value-of select="substring-before(attributes/attribute[@key='transcript_id']/text(),'.')"/>
			</xsl:attribute>
		</xsl:when>
		<xsl:otherwise><!-- anonymous no rdf:about --></xsl:otherwise>
	</xsl:choose>

	<xsl:apply-templates select="attributes"/>

	<!-- add coordinates -->
	<xsl:apply-templates select="." mode="coord"/>

</xsl:element>

</xsl:template>



<!--===========================================================================-->


<!--===========================================================================-->

<xsl:template match="attributes" >
<xsl:apply-templates select="attribute"/>
</xsl:template>


<xsl:template match="attribute" >
	<xsl:choose>
		<xsl:when test="(@key='gene_biotype' or @key='gene_type') and ../../name() = 'gene'">
			<xsl:apply-templates select="." mode="biotype"/>
		</xsl:when>
		<xsl:when test="@key='transcript_biotype' and ../../name() = 'transcript'">
			<xsl:apply-templates select="." mode="biotype"/>
		</xsl:when>
		<xsl:when test="@key='gene_id'">
			<xsl:if test="../../name() != 'gene'">
				<u:has_gene>
					<xsl:attribute name="u:has_gene">
						<xsl:text>http://rdf.ebi.ac.uk/resource/ensembl/</xsl:text>
						<xsl:value-of select="substring-before(./text(),'.')"/>
					</xsl:attribute>
				</u:has_gene>
			</xsl:if>
		</xsl:when>
		<xsl:when test="@key='transcript_id'">
			<xsl:if test="../../name() != 'transcript'">
					<u:has_gene>
						<xsl:attribute name="u:has_transcript">
							<xsl:text>http://rdf.ebi.ac.uk/resource/ensembl/</xsl:text>
							<xsl:value-of select="substring-before(./text(),'.')"/>
						</xsl:attribute>
					</u:has_gene>
				</xsl:if>
		</xsl:when>
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
<xsl:template match="attribute" mode="biotype" >


<xsl:variable name="typename">
<xsl:call-template name="biotype2qname">
  <xsl:with-param name="type" select="text()"/>
  <xsl:with-param name="default" ><xsl:text></xsl:text></xsl:with-param>
</xsl:call-template>
</xsl:variable>

<xsl:if test="$typename !="" ">
<rdf:type>
	<xsl:attribute name="rdf:resource">
		<xsl:text>&obo;</xsl:text>
		<xsl:value-of select="$typename"/>
	</xsl:attribute>
</rdf:type>
</xsl:if>



</xsl:template>
<!--===========================================================================-->


<xsl:template match="*" mode="coord" >
<u:has_interval>
<u:Interval>
<u:chromosome>
<xsl:attribute name="rdf:resource">
	<xsl:text>&u;build/</xsl:text>
	<xsl:value-of select="/gtf/dictionary/@build"/>
	<xsl:text>/</xsl:text>
	<xsl:value-of select="@chrom"/>
</xsl:attribute>
</u:chromosome>
<u:start rdf:datatype="&xsd;int"><xsl:value-of select="@start"/></u:start>
<u:end rdf:datatype="&xsd;int"><xsl:value-of select="@end"/></u:end>
</u:Interval>
</u:has_interval>
</xsl:template>


</xsl:stylesheet>
