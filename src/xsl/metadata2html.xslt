<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet
        version="1.0"
        xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
 >
<xsl:output method="xml" encoding="UTF-8"  indent="yes" omit-xml-declaration="yes"/>

<!-- ========================================================================================== -->
<!-- ========================================================================================== -->
<!-- ========================================================================================== -->

<xsl:template match="/">

<xsl:apply-templates select="schema"/>

</xsl:template>

<!-- ========================================================================================== -->
<!-- ========================================================================================== -->
<!-- ========================================================================================== -->

<xsl:template match="schema">
<html>
<head>
  <title><xsl:apply-templates select="@name"/></title>
</head>
<body>
<h1><xsl:apply-templates select="@name"/></h1>
<xsl:apply-templates select="description"/>
<xsl:apply-templates select="params"/>

<h2>Custom Tool Arguments</h2>

<p>A pipeline might not always support every possible argument or option of a particular tool used in pipeline.
Fortunately, nextflow  some freedom to users to insert additional parameters that the pipeline does not include by default.
</p>

<h2>Credits</h2>
<div>
</div>
<h2>Contributions and Support</h2>
<div>
If you would like to contribute to this pipeline, please submit a pull-request
</div>

</body>
</html>
</xsl:template>
<!-- ========================================================================================== -->
<!-- ========================================================================================== -->
<!-- ========================================================================================== -->
<xsl:template match="description">
<p class="desc">
<xsl:apply-templates select="*|text()"/>
</p>
</xsl:template>
<!-- ========================================================================================== -->
<!-- ========================================================================================== -->
<!-- ========================================================================================== -->


<xsl:template match="params">
<h2>Parameters</h2>
<xsl:apply-templates select="section"/>
</xsl:template>

<!-- ========================================================================================== -->
<!-- ========================================================================================== -->
<!-- ========================================================================================== -->


<xsl:template match="section">
<h3><xsl:apply-templates select="@name"/></h3>
<xsl:apply-templates select="description"/>
<xsl:apply-templates select="help_text"/>

<table>
<thead>
  <caption><xsl:apply-templates select="@name"/></caption>
  <tr>
    <th>name</th>
    <th>description</th>
    <th>type</th>
    <th>default</th>
    <th>format</th>
    <th>pattern</th>
    <th>help</th>
  </tr>
</thead>
<tbody>
  <xsl:apply-templates select="param"/>
</tbody>
</table>
</xsl:template>

<!-- ========================================================================================== -->
<!-- ========================================================================================== -->
<!-- ========================================================================================== -->

<xsl:template match="param">
<tr>
  <td><code>--<xsl:apply-templates select="@name"/></code></td>
  <td><xsl:apply-templates select="description"/></td>
  <td><xsl:apply-templates select="type"/></td>
  <td><xsl:apply-templates select="default"/></td>
  <td><xsl:apply-templates select="format"/></td>
  <td><xsl:apply-templates select="pattern"/></td>
  <td><xsl:apply-templates select="help_text"/></td>
</tr>
</xsl:template>



<!-- ========================================================================================== -->
<!-- ========================================================================================== -->
<!-- ========================================================================================== -->

</xsl:stylesheet>
