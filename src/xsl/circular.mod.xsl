<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet
 version="1.0"
 xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
 xmlns:h="http://www.w3.org/1999/xhtml"
 xmlns="http://www.w3.org/2000/svg"
 xmlns:math="http://exslt.org/math"
 extension-element-prefixes="math"
 >



<xsl:variable name="PI" select="number(3.14159265359)"/>
<xsl:variable name="lowercase" select="'abcdefghijklmnopqrstuvwxyz'" />
<xsl:variable name="uppercase" select="'ABCDEFGHIJKLMNOPQRSTUVWXYZ'" />

<!-- ==========================================================-->

<xsl:template name="donut">
<xsl:param name="outer_radius"/>
<xsl:param name="inner_radius"/>
<xsl:param name="class"></xsl:param>
<xsl:param name="style"></xsl:param>
<xsl:param name="title"></xsl:param>

<path  fill-rule="evenodd">
 <xsl:attribute name="d">
 
	<xsl:call-template name="donut_path">
			<xsl:with-param name="outer_radius" select="$outer_radius"/>
			<xsl:with-param name="inner_radius" select="$inner_radius"/>
	</xsl:call-template>

 </xsl:attribute>

	<xsl:if test="string-length(normalize-space($class)) &gt; 0">
		 <xsl:attribute name="class">
		 	<xsl:value-of select="$class"/>
		 </xsl:attribute>
	</xsl:if>

	<xsl:if test="string-length(normalize-space($style)) &gt; 0">
		 <xsl:attribute name="style">
		 	<xsl:value-of select="$style"/>
		 </xsl:attribute>
	</xsl:if>

	<xsl:if test="string-length(normalize-space($title)) &gt; 0">
		<title><xsl:value-of select="$title"/></title>
	</xsl:if>
</path>

</xsl:template>



<!-- ==========================================================-->

<xsl:template name="donut_path">
<xsl:param name="outer_radius"/>
<xsl:param name="inner_radius"/>

	<xsl:call-template name="circle_as_path">
		<xsl:with-param name="r" select="$outer_radius"/>
	</xsl:call-template>
	<xsl:call-template name="circle_as_path">
		<xsl:with-param name="r" select="$inner_radius"/>
	</xsl:call-template>

</xsl:template>


<!-- ==========================================================-->
<xsl:template name="circle_as_path">
<xsl:param name="r"/>
	<xsl:text> M </xsl:text>
	<xsl:value-of select="$r"/>

	<xsl:text>  0 A </xsl:text>
	<xsl:value-of select="$r"/>
	<xsl:text>  </xsl:text>
	<xsl:value-of select="$r"/>
	<xsl:text>  0 1 0 </xsl:text>
	<xsl:value-of select="number($r) * -1"/>
	<xsl:text>  0 A  </xsl:text>

	<xsl:value-of select="$r"/>
	<xsl:text>  </xsl:text>
	<xsl:value-of select="$r"/>
	<xsl:text>  0 1 0 </xsl:text>
	<xsl:value-of select="number($r)"/>
	<xsl:text>  0 Z </xsl:text>
</xsl:template>


<!-- ==========================================================-->
<xsl:template name="arc">
<xsl:param name="f1"/>
<xsl:param name="f2"/>
<xsl:param name="outer_radius"/>
<xsl:param name="inner_radius"/>

<xsl:text>M </xsl:text>

<xsl:call-template name="toPoint">
	<xsl:with-param name="f">
		<xsl:value-of select="$f1"/>
	</xsl:with-param>
	<xsl:with-param name="r">
		<xsl:value-of select="$inner_radius"/>
	</xsl:with-param>
</xsl:call-template>

<xsl:text> L </xsl:text>

<xsl:call-template name="toPoint">
	<xsl:with-param name="f">
		<xsl:value-of select="$f1"/>
	</xsl:with-param>
	<xsl:with-param name="r">
		<xsl:value-of select="$outer_radius"/>
	</xsl:with-param>
</xsl:call-template>

<xsl:call-template name="arcTo">
	<xsl:with-param name="r">
		<xsl:value-of select="$outer_radius"/>
	</xsl:with-param>
	<xsl:with-param name="f1">
		<xsl:value-of select="$f1"/>
	</xsl:with-param>
	<xsl:with-param name="f2">
		<xsl:value-of select="$f2"/>
	</xsl:with-param>
</xsl:call-template>

<xsl:text> L </xsl:text>

<xsl:call-template name="toPoint">
	<xsl:with-param name="f">
		<xsl:value-of select="$f2"/>
	</xsl:with-param>
	<xsl:with-param name="r">
		<xsl:value-of select="$inner_radius"/>
	</xsl:with-param>
</xsl:call-template>

<xsl:call-template name="arcTo">
	<xsl:with-param name="r">
		<xsl:value-of select="$inner_radius"/>
	</xsl:with-param>
	<xsl:with-param name="f1">
		<xsl:value-of select="$f2"/>
	</xsl:with-param>
	<xsl:with-param name="f2">
		<xsl:value-of select="$f1"/>
	</xsl:with-param>
</xsl:call-template>

 

</xsl:template>

<!-- convert fraction/radius to X(space)Y coordinate -->
<xsl:template name="toPoint">
<xsl:param name="f"/>
<xsl:param name="r"/>
<xsl:variable name="angle" select="(($PI * 2) * number($f)) - ($PI div 2)" />
<xsl:value-of select="number($r)*math:cos($angle)"/>
<xsl:text> </xsl:text>
<xsl:value-of select="number($r)*math:sin($angle)"/>
</xsl:template>

<!-- convert fraction/radius to X coordinate -->
<xsl:template name="toX">
<xsl:param name="f"/>
<xsl:param name="r"/>
<xsl:variable name="angle" select="(($PI * 2) * number($f)) - ($PI div 2)" />
<xsl:value-of select="number($r)*math:cos($angle)"/>
</xsl:template>

<!-- convert fraction/radius to Y coordinate -->
<xsl:template name="toY">
<xsl:param name="f"/>
<xsl:param name="r"/>
<xsl:variable name="angle" select="(($PI * 2) * number($f)) - ($PI div 2)" />
<xsl:value-of select="number($r)*math:sin($angle)"/>
</xsl:template>

<xsl:template name="arcTo">
<xsl:param name="r"/>
<xsl:param name="f1"/>
<xsl:param name="f2"/>

<xsl:text> A </xsl:text>
<xsl:value-of select="$r"/>
<xsl:text> </xsl:text>
<xsl:value-of select="$r"/>
<xsl:text> 0 0 </xsl:text>
<xsl:choose>
	<xsl:when test="number($f1) &lt; number($f2) ">
		<xsl:text>1 </xsl:text>
	</xsl:when>
	<xsl:otherwise>
		<xsl:text>0 </xsl:text>
	</xsl:otherwise>
</xsl:choose>

<xsl:call-template name="toPoint">
	<xsl:with-param name="f">
		<xsl:value-of select="$f2"/>
	</xsl:with-param>
	<xsl:with-param name="r">
		<xsl:value-of select="$r"/>
	</xsl:with-param>
</xsl:call-template>

</xsl:template>

<xsl:template name="x">
<xsl:variable name="angle" select="(($PI * 2) * number($f)) - ($PI div 2)" />
<xsl:value-of select="number($r)*math:cos($angle)"/>
<xsl:text> </xsl:text>
<xsl:value-of select="number($r)*math:sin($angle)"/>
</xsl:template>


<xsl:template name="cytoband_color">
<xsl:param name="stain"/>
<xsl:choose>
	<xsl:when test='starts-with($stain,"gneg")'>
		<xsl:text>lightblue</xsl:text>
	</xsl:when>
	<xsl:when test='starts-with($stain,"acen")'>
		<xsl:text>orange</xsl:text>
	</xsl:when>
	<xsl:when test='starts-with($stain,"gvar")'>
		<xsl:text>slategray</xsl:text>
	</xsl:when>
	<xsl:when test='starts-with($stain,"gpos") and string-length($stain) &gt; 4'>
		<xsl:variable name="percent1" select="number(substring-after($stain,'s'))" />
		<xsl:variable name="gray" select="format-number( 50 + (215 * ($percent1 div 100.0)) , '#')"/>
		<xsl:text>rgb(</xsl:text>
		<xsl:value-of select="$gray"/>
		<xsl:text>,</xsl:text>
		<xsl:value-of select="$gray"/>
		<xsl:text>,</xsl:text>
		<xsl:value-of select="$gray"/>
		<xsl:text>)</xsl:text>
	</xsl:when>
	<xsl:otherwise>
		<xsl:text>honeydew</xsl:text>
	</xsl:otherwise>
    </xsl:choose>
</xsl:template>


<xsl:template name="ucsc_anchor">
<xsl:param name="build"/>
<xsl:param name="contig"/>
<xsl:param name="start"/>
<xsl:param name="end"/>
<xsl:variable name="b2" select="translate($build, $lowercase, $uppercase)"/>
<xsl:variable name="n">
	<xsl:choose>
		<xsl:when test="$b2='HG19' or $b2='GRCH37' ">
			<xsl:text>org=Human&amp;db=hg19</xsl:text>
		</xsl:when>
		<xsl:when test="$b2='HG38' or $b2='GRCH38' ">
			<xsl:text>org=Human&amp;db=hg38</xsl:text>
		</xsl:when>
		<xsl:otherwise>
			<xsl:text></xsl:text>
		</xsl:otherwise>
	</xsl:choose>
</xsl:variable>

<xsl:choose>
	<xsl:when test="string-length($n) &gt; 0">
		<xsl:text>http://genome.ucsc.edu/cgi-bin/hgTracks?</xsl:text>
		<xsl:value-of select="$n"/>
		<xsl:text>&amp;position=</xsl:text>
		<xsl:value-of select="$contig"/>
		<xsl:text>%3A</xsl:text>
		<xsl:value-of select="$start"/>
		<xsl:text>-</xsl:text>
		<xsl:value-of select="$end"/>
	</xsl:when>
	<xsl:otherwise>
		<xsl:text>#</xsl:text>
	</xsl:otherwise>
</xsl:choose>

</xsl:template>

</xsl:stylesheet>
