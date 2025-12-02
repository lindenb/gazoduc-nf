<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet
 version="1.0"
 xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
 xmlns="http://www.w3.org/2000/svg"
 >
<!--

this stylesheet converts a faix+jvarkit bed2xml to a SVG fragment

-->


<xsl:output method="xml" encoding="UTF-8" />
<xsl:include href="circular.mod.xsl"/>
<xsl:param name="radius">1000</xsl:param>

<xsl:template match="/">
<xsl:apply-templates select="bed/body"/>
</xsl:template>

<xsl:template match="body">
<svg width="2100" height="2100">
<style>
path.contig {
		fill: gray;
		fill-opacity:0.3;
		}
path.contigY {
		fill: aliceblue;
		}
path.contigX {
		fill: lightpink;
		}
path.contigmod1 {
		fill:lavender;
		}
path.contigmod0 {
		fill:lavenderblush;
		}
line.contig {
		stroke: gray;
		stroke-width:1px;
		stroke-dasharray:2px, 5px;
		}
text.contig {
		text-anchor:middle;
		}
</style>
<g transform="translate(1000,1000)">
<xsl:apply-templates select="rec"/>
</g>
</svg>
</xsl:template>


<xsl:template match="rec">


<xsl:variable name="r1" select="number($radius)"/>
<xsl:variable name="r2" select="number($radius) * 0.001"/>

<g>
 <xsl:attribute name="id">
 	 <xsl:text>contig.</xsl:text>
 	 <xsl:value-of select="@contig"/>
 </xsl:attribute>
<path>
     <xsl:attribute name="class">
        <xsl:choose>
                <xsl:when test='@contig="chrX" or @contig="X"'>
                        <xsl:text>contig contigX</xsl:text>
                </xsl:when>
                <xsl:when test='@contig="chrY" or @contig="Y"'>
                        <xsl:text>contig contigY</xsl:text>
                </xsl:when>
                <xsl:when test='position() mod 2 = 0'>
                        <xsl:text>contig contigmod0</xsl:text>
                </xsl:when>
                <xsl:otherwise>
                        <xsl:text>contig contigmod1</xsl:text>
                </xsl:otherwise>
        </xsl:choose>
    </xsl:attribute>
    
    
    <xsl:attribute name="d">
        <xsl:call-template name="arc">
                <xsl:with-param name="f1">
                        <xsl:value-of select="number(@f1)"/>
                </xsl:with-param>
                <xsl:with-param name="f2">
                        <xsl:value-of select="number(@f2)"/>
                </xsl:with-param>
                <xsl:with-param name="outer_radius">
                        <xsl:value-of select="$r1"/>
                </xsl:with-param>
                <xsl:with-param name="inner_radius">
                        <xsl:value-of select="$r2"/>
                </xsl:with-param>
        </xsl:call-template>
	</xsl:attribute>
	<title>
		<xsl:value-of select="@contig"/>
		<xsl:text> </xsl:text>
		<xsl:value-of select="@length"/>
		<xsl:text>bp</xsl:text>
	</title>
</path>
	
	<!-- left line to center -->
	<line class="contig">
		<!-- x1 -->
		<xsl:attribute name="x1">
		 <xsl:call-template name="toX">
                	<xsl:with-param name="f">
                        	<xsl:value-of select="number(@f1)"/>
                	</xsl:with-param>
                	<xsl:with-param name="r">
                        	<xsl:value-of select="$r1 + 10.0"/>
                	</xsl:with-param>
        	</xsl:call-template>
		</xsl:attribute>
		<!-- y1 -->
		<xsl:attribute name="y1">
		 <xsl:call-template name="toY">
                	<xsl:with-param name="f">
                        	<xsl:value-of select="number(@f1)"/>
                	</xsl:with-param>
                	<xsl:with-param name="r">
                        	<xsl:value-of select="$r1 + 10.0"/>
                	</xsl:with-param>
        	</xsl:call-template>
		</xsl:attribute>
		<!-- x2 -->
		<xsl:attribute name="x2">
		 <xsl:call-template name="toX">
                	<xsl:with-param name="f">
                        	<xsl:value-of select="number(@f1)"/>
                	</xsl:with-param>
                	<xsl:with-param name="r">
                        	<xsl:value-of select="$r2"/>
                	</xsl:with-param>
        	</xsl:call-template>
		</xsl:attribute>
		<!-- y2 -->
		<xsl:attribute name="y2">
		 <xsl:call-template name="toY">
                	<xsl:with-param name="f">
                        	<xsl:value-of select="number(@f1)"/>
                	</xsl:with-param>
                	<xsl:with-param name="r">
                        	<xsl:value-of select="$r2"/>
                	</xsl:with-param>
        	</xsl:call-template>
		</xsl:attribute>
	</line>
   
   <!-- right line to center -->
  
	<line class="contig">
		<!-- x1 -->
		<xsl:attribute name="x1">
		 <xsl:call-template name="toX">
                	<xsl:with-param name="f">
                        	<xsl:value-of select="number(@f2)"/>
                	</xsl:with-param>
                	<xsl:with-param name="r">
                        	<xsl:value-of select="$r1 + 30.0"/>
                	</xsl:with-param>
        	</xsl:call-template>
		</xsl:attribute>
		<!-- y1 -->
		<xsl:attribute name="y1">
		 <xsl:call-template name="toY">
                	<xsl:with-param name="f">
                        	<xsl:value-of select="number(@f2)"/>
                	</xsl:with-param>
                	<xsl:with-param name="r">
                        	<xsl:value-of select="$r1 + 30.0"/>
                	</xsl:with-param>
        	</xsl:call-template>
		</xsl:attribute>
		<!-- x2 -->
		<xsl:attribute name="x2">
		 <xsl:call-template name="toX">
                	<xsl:with-param name="f">
                        	<xsl:value-of select="number(@f2)"/>
                	</xsl:with-param>
                	<xsl:with-param name="r">
                        	<xsl:value-of select="$r2"/>
                	</xsl:with-param>
        	</xsl:call-template>
		</xsl:attribute>
		<!-- y2 -->
		<xsl:attribute name="y2">
		 <xsl:call-template name="toY">
                	<xsl:with-param name="f">
                        	<xsl:value-of select="number(@f2)"/>
                	</xsl:with-param>
                	<xsl:with-param name="r">
                        	<xsl:value-of select="$r2"/>
                	</xsl:with-param>
        	</xsl:call-template>
		</xsl:attribute>
	</line>  
  
  
	<text class="contig">
		<xsl:attribute name="x">0</xsl:attribute>
		<xsl:attribute name="y"><xsl:value-of select="$r1 + 30"/></xsl:attribute>
		<xsl:attribute name="transform">
			<xsl:text>rotate(</xsl:text>
			<!-- +180 works but I don't understand why ... -->
			<xsl:value-of select="((number(@f1) + number(@f2) ) div 2.0) * 360.0 + 180"/>
        	<xsl:text>)</xsl:text>
		</xsl:attribute>
		
		
		<xsl:value-of select="@contig"/>
	</text>

</g>
</xsl:template>

</xsl:stylesheet>
