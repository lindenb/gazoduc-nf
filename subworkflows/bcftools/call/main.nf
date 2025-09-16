/*

Copyright (c) 2024 Pierre Lindenbaum

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
include {BCFTOOLS_CALL as CALL } from '../../../modules/bcftools/call'
include {BCFTOOLS_MERGE        } from '../../../modules/bcftools/merge'
include {BCFTOOLS_CONCAT       } from '../../../modules/bcftools/concat'



workflow BCFTOOLS_CALL {
	take:
		meta
		fasta
		fai
		dict
		pedigree
		beds
		bams //[ meta, bam,bai]
	main:
		versions = Channel.empty()
		ch1 = bams
			.map{
				if(it[0].containsKey("batch")) return it;
				return [it[0].plus(batch:it[0].sample),it[1],it[2]];
				}
			.map{[[id:it[0].batch],[it[1],it[2]]]}
			.groupTuple()
			.map{[it[0],it[1].flatten()]}
		
		ch2 = ch1.combine(beds.map{it[1]})

		CALL(
			fasta,
			fai,
			pedigree,
			[[id:"noploidy"],[]],
			ch2
			)
		versions = versions.mix(CALL.out.versions)

		ch3 = CALL.out.vcf
				.map{[[id:it[3].toRealPath().toString()/* bed */],[it[1],it[2]]]}
				.groupTuple()
				.map{[it[0],it[1].flatten(), [] /* no bed */]}
				.branch{v->
					to_merge:v[1].size()>2 /* bed and it's index */
					other: true
					}
		// merge thos having more than one sample
		BCFTOOLS_MERGE( ch3.to_merge )
		versions = versions.mix(BCFTOOLS_MERGE.out.versions)

		without_merge = ch3.other.map{[
			it[0],
			it[1].find{v->v.name.endsWith(".bcf") || v.name.endsWith("*.vcf.gz")}, 
			it[1].find{v->v.name.endsWith(".tbi") || v.name.endsWith("*.csi")}
			]}

		SET_GQ(BCFTOOLS_MERGE.out.vcf.mix(without_merge))
		versions = versions.mix(SET_GQ.out.versions)


		BCFTOOLS_CONCAT(
			SET_GQ.out.vcf
				.map{[[id:"call"],[it[1],it[2]]]}
				.groupTuple()
				.map{[it[0],it[1].flatten()]},
			[[id:"nobed"],[]]
			)
		versions = versions.mix(BCFTOOLS_CONCAT.out.versions)

	emit:
		versions
		vcf = BCFTOOLS_CONCAT.out.vcf
	}




process SET_GQ {
tag "${meta.id?:""}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
   tuple val(meta),path(vcf),path(vcfidx)
output:
    tuple val(meta),path("*.bcf"),path("*.csi"),emit:vcf
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"\${MD5}"+".gq"
"""
set -x
mkdir -p TMP
MD5=`cat "${vcf}" | md5sum | cut -d ' ' -f1`

cat << __EOF__ > TMP/jeter.code
final VariantContextBuilder vcb = new VariantContextBuilder(variant);
vcb.genotypes(
  variant.getGenotypes().stream().map(G->{
	if(G.hasGQ()) return G;
        final GenotypeBuilder gb=new GenotypeBuilder(G);
        final int dp= (G.hasDP()?G.getDP():30);
        double gq = 99;
        if(dp<25) gq = gq * (dp/25.0);
        if(G.isNoCall()) {
                gq=0;
                }
        else if(G.hasAD() && G.getAD().length==2) {
                final int[] ad=G.getAD();
                double dp2 = ad[0]+ad[1];
                if((G.isHomRef() || G.isHomVar()) && dp2>0) {
                        gq = gq * Math.max(ad[0],ad[1])/dp2;
                        }
                else if(G.isHet() && dp2>0) {
                        gq = gq * ((Math.min(ad[0],ad[1])/dp2) / 0.5 );
                        }
                }
        gb.GQ((int)gq);
        return gb.make();
        }).collect(Collectors.toList())
  );
return vcb.make();
__EOF__


bcftools view "${vcf}" |\\
	awk '/^#CHROM/ {printf("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\\"Genotype Quality for bcftools\\">\\n");} {print;}' |\\
	jvarkit  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP vcffilterjdk -f ~/jeter.code |\\
	bcftools view -O b -o TMP/jeter.bcf

bcftools index --threads ${task.cpus} -f TMP/jeter.bcf

mv TMP/jeter.bcf     ${prefix}.bcf
mv TMP/jeter.bcf.csi ${prefix}.bcf.csi

cat << END_VERSIONS > versions.yml
${task.process}:
    bcftools: \$(bcftools version | awk '(NR==1)  {print \$NF}')
    jvarkit: todo
END_VERSIONS
"""
}
