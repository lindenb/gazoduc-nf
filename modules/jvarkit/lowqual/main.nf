/*

Copyright (c) 2026 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
The MIT License (MIT)
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

process JVARKIT_FILTER_LOWQUAL {
	label "process_single"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	afterScript "rm -rf TMP"
	tag "${meta.id?:vcf.name}"
	input:
		tuple val(meta ),path(vcf)
	output:
		tuple val(meta ),path("*.vcf.gz"),emit:vcf
		path("versions.yml"),emit:versions
	script:
		def args1 = task.ext.args1?:""
		def args2 = task.ext.args2?:""
		def prefix  = task.ext.prefix?:"${meta.id}.lowqual"
        
		def minRatioSingleton = task.ext.ad_ratio?:0.2
        def minGQ = task.ext.minGQ?:60
        def minDP = task.ext.minDP?:10
        def maxDP = task.ext.maxDP?:300
	"""
	hostname 1>&2
	set -o pipefail

	mkdir -p TMP

cat << EOF > TMP/jeter.code
final  VariantContextBuilder vcb = new VariantContextBuilder(variant);

if(  variant.getGenotypes().stream().anyMatch(G->G.isHet() && G.hasAD()) &&
     variant.getGenotypes().stream()
    .filter(G->G.isHet() && G.hasAD())
    .map(G->G.getAD())
    .filter(A->A.length>=2)
    .allMatch(A0->{
	// if more than 2 item in AD heterozgous, sort and keep the two highest
	final int[] A= Arrays.copyOf(A0, A0.length);
	Arrays.sort(A);
	final int a0 = A[A.length-2];
	final int a1 = A[A.length-1];
        final double n = a0 + a1;
        if(n<=0) return true;
        double r= a1/n;
        if(r < ${minRatioSingleton} || r> (1.0 - ${minRatioSingleton})) return true;
        return false;
        })) {
        vcb.filter("HET_BAD_AD_RATIO");
        }

if(variant.getGenotypes().stream()
    .filter(G->G.hasAltAllele() && G.hasDP())
    .mapToInt(G->G.getDP())
    .allMatch(DP->DP < ${minDP})) {
        vcb.filter("LOW_DEPTH");
        }

if(variant.getGenotypes().stream()
    .filter(G->G.hasAltAllele() && G.hasDP())
    .mapToInt(G->G.getDP())
    .allMatch(DP->DP > ${maxDP})) {
        vcb.filter("HIGH_DEPTH");
        }

if(variant.getGenotypes().stream()
    .filter(G->G.hasAltAllele() && G.hasGQ())
    .mapToInt(G->G.getGQ())
    .allMatch(GQ->GQ < ${minGQ})) {
        vcb.filter("LOW_GQ");
        }

if(variant.hasAttribute("FS") && variant.getAttributeAsDouble("FS",0) > 60 ) {
    vcb.filter("HIGH_FS");
    }


if(variant.hasAttribute("QD") && variant.getAttributeAsDouble("QD",10) < 2.0 ) {
    vcb.filter("LOW_QD");
    }


if(variant.hasAttribute("SOR") && variant.getAttributeAsDouble("SOR",0) > 3.0 ) {
    vcb.filter("HIGH_SOR");
    }

if(variant.hasAttribute("ReadPosRankSum") && variant.getAttributeAsDouble("ReadPosRankSum",0) < -8.0 ) {
    vcb.filter("LOW_ReadPosRankSum");
    }

if(variant.hasAttribute("MQ") && variant.getAttributeAsDouble("MQ",100) < 40.0 ) {
    vcb.filter("LOW_MQ");
    }

if(variant.hasAttribute("MQRankSum") && variant.getAttributeAsDouble("MQRankSum",0) < -12.5 ) {
    vcb.filter("LOW_MQRankSum");
    }

return vcb.make();
EOF

bcftools view  ${args1} "${vcf}" |\\
	jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP vcffilterjdk \\
            --extra-filters 'LOW_MQRankSum,LOW_MQ,LOW_ReadPosRankSum,HIGH_FS,HIGH_SOR,LOW_QD,LOW_GQ,LOW_DEPTH,HIGH_DEPTH,HET_BAD_AD_RATIO' \\
            -f TMP/jeter.code |\\
	bcftools view ${args2}  -O z -o TMP/${prefix}.vcf.gz

	mv TMP/${prefix}.vcf.gz ./

cat << END_VERSIONS > versions.yml
"${task.process}":
	jvarkit: todo
END_VERSIONS
	"""
stub:
def prefix  = task.ext.prefix?:"${meta.id}.lowqual"
"""
touch versions.yml ${prefix}.vcf.gz
"""
}
