/*

Copyright (c) 2025 Pierre Lindenbaum

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
include {SMOOVE_CALL                        } from '../../modules/smoove/call'
include {SMOOVE_MERGE                       } from '../../modules/smoove/merge'
include {SMOOVE_GENOTYPE                    } from '../../modules/smoove/genotype'
include {SMOOVE_PASTE                       } from '../../modules/smoove/paste'

workflow SMOOVE {
	take:
		meta
		fasta
		fai
		dict
		exclude_bed
		bams
	main:

		versions = Channel.empty()

        /* if there is no meta.status, treat everyone as case */
        bams.map{[it[0].status ? it[0] : it[0].plus("status":"case"), it[1], it[2]] }.
			branch{
				controls : it[0].status && it[0].status.equals("control")
				cases : true
			}.set{bams_status_ch}


		SMOOVE_CALL(
			fasta,
			fai,
			exclude_bed,
			bams_status_ch.cases
			)
		versions = versions.mix(SMOOVE_CALL.out.versions)
	

		ch1 = SMOOVE_CALL.out.vcf
				.map{[it[1],it[2]]}
				.collect()
				.map{[meta,it]}

		
		SMOOVE_MERGE(
			fasta,
			fai,
			ch1
			)

		versions = versions.mix(SMOOVE_MERGE.out.versions)	
		

		SMOOVE_GENOTYPE(
			fasta,
			fai,
			SMOOVE_MERGE.out.vcf,
			bams
		)

		versions = versions.mix(SMOOVE_GENOTYPE.out.versions)

		SMOOVE_PASTE(fasta,fai,dict,SMOOVE_GENOTYPE.out.vcf
			.map{[it[1],it[2]]}
			.collect()
			.map{[meta,it]}
		)
		versions = versions.mix(SMOOVE_PASTE.out.versions)
	
	emit:
		versions
		vcf = SMOOVE_PASTE.out.vcf
	}




process PASTE_ALL {
tag "N=${L.size()}"
cache "lenient"
afterScript  "rm -rf TMP TMP2"
memory "10g"
input:
	val(meta)
	val(reference)
	val(img)
	val(L)
output:
	path("${params.prefix?:""}smoove.bcf"),emit:vcf
	path("${params.prefix?:""}smoove.bcf.csi"),emit:index
	path("version.xml"),emit:version
script:
	def prefix = params.prefix?:""
	def gff = ""//gff3
	log.warn("JE SUPPRIME GFF POUR LE MOMENT. CA BUG POUR SOLENA OCT 2022")
"""
	hostname 1>&2
	${moduleLoad("bcftools picard")}
	set -x
	mkdir -p TMP TMP2
	# smoove will write to the system TMPDIR. For large cohorts, make sure to set this to something with a lot of space
	export TMPDIR=\${PWD}/TMP2

if [ "${L.size()}" -ne "1" ] ; then

cat << EOF > TMP/jeter.list
${L.join("\n")}
EOF

	#module load singularity/2.4.5
	singularity exec\
		--home \${PWD} \
		--bind \${PWD}/TMP:/outdir \
			`xargs -a TMP/jeter.list -L 1 dirname | awk '{printf(" --bind %s:/d%d ",\$0,NR);}'  ` \
		${img} \
		smoove paste --name paste  \
			 `awk -F '/' '{printf(" /d%d/%s ",NR,\$NF);}' TMP/jeter.list`

	mv paste.smoove.square.vcf.gz TMP/jeter1.vcf.gz

else

	bcftools view -O z -o TMP/jeter1.vcf.gz "${L[0]}"

fi


	if [ ! -z "${gff}" ] ; then
	
		singularity exec\
			--home \${PWD} \
			--bind \${PWD}/TMP:/outdir \
			${gff.isEmpty()?"":" --bind "+ file(gff).getParent()+":/gdir"} \
			--bind \${PWD}/TMP:/data1 \
			${img} \
			smoove annotate ${gff.isEmpty()?"":"--gff /gdir/"+file(gff).name} /data1/jeter1.vcf.gz > TMP/jeter2.vcf

		bcftools view -O z -o TMP/jeter2.vcf.gz TMP/jeter2.vcf

		mv TMP/jeter2.vcf.gz TMP/jeter1.vcf.gz
	fi


	#######################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">smoove paste all and annotate</entry>
		<entry key="gff">${gff}</entry>
		<entry key="count">${L.size()}</entry>
		<entry key="versions">${getVersionCmd("bcftools picard/UpdateVcfSequenceDictionary")}</entry>
	</properties>
	EOF
"""
}
