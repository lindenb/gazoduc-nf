include {moduleLoad} from '../utils/functions.nf'

process RVTESTS01_VCF_01 {
tag "${file(vcf).name}"
afterScript "rm -rf TMP"
memory "10g"
input:
	val(meta)
	val(reference)
	val(vcf)
	val(pedigree)
output:
	path("assoc.list"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("rvtests bcftools jvarkit")}

mkdir TMP ASSOC

## https://github.com/zhanxw/rvtests/issues/80
cut -f 1 "${reference}.fai" > TMP/chroms.A.txt
sed 's/^chr//' TMP/chroms.A.txt > TMP/chroms.B.txt
paste TMP/chroms.A.txt TMP/chroms.B.txt > TMP/chroms.C.txt


bcftools annotate -O v --rename-chrs TMP/chroms.C.txt "${vcf}" |\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/vcfgenesplitter.jar \
		-E 'ANN/GeneId ANN/FeatureId VEP/GeneId VEP/Ensp VEP/Feature' \
		--manifest TMP/jeter.mf \
		-o \${PWD}/TMP

i=1
	awk -F '\t' '/^[^#]/ {S=\$4;gsub(/[\\/]+/,"_",S);G=\$5;gsub(/[\\/_\\-]+/,"_",G);K=\$6;gsub(/[\\/]+/,"_",K);printf("%s_%s_%s\t%s:%d-%d\t%s\\n",K,G,S,\$1,\$2,\$3,\$7);}' TMP/jeter.mf |\
	while read K RANGE V
	do
		bcftools index -t "\${V}"

		# build setFile
		echo "\$K\t\${RANGE}" > TMP/variants.setfile

		rvtest  --noweb \
        		--inVcf "\${V}" \
			--setFile TMP/variants.setfile \
	        	--pheno "${pedigree}" \
		        --out "ASSOC/part.\${i}" \
        		--burden cmc,zeggini,mb,fp,exactCMC,cmcWald,rarecover,cmat \
	        	--vt price,analytic \
		        --kernel 'skat[nPerm=1000],kbac,skato'

		i=\$((i+1))
	done


find \${PWD}/ASSOC -type f -name "part.*assoc" > assoc.list

cat << EOF > version.xml
<properties id="${task.process}">
  <entry key="name">${task.process}</entry>
  <entry key="description">VCF is split by transcript / gene, and then rvtest is applied.</entry>
  <entry key="name"Pedigree">${pedigree}</entry>
</properties>
EOF
"""
stub:
"""
touch assoc.list
echo "<properties/>" > version.xml
"""
}
