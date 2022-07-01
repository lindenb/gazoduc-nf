
process PEDIGREE_FOR_RVTESTS {
executor "local"
input:
        val(meta)
        val(cases)
        val(controls)
output:
        path("rvtests.pedigree.ped"),emit:pedigree
        path("version.xml"),emit:version
script:
	def p=meta.rvtests_phenotype_name?:"y1"
"""
set -o pipefail

# assert no duplicate between cases and controls
cat "${cases}" "${controls}" | sort | uniq -d  > dups.txt
test ! -s dups.txt
rm dups.txt

echo "fid\tiid\tfatid\tmatid\tsex\t${p}" > jeter.ped
awk '{printf("%s\t%s\t0\t0\t0\t2\\n",\$1,\$1);}' "${cases}"    >> jeter.ped
awk '{printf("%s\t%s\t0\t0\t0\t1\\n",\$1,\$1);}' "${controls}" >> jeter.ped

mv -v jeter.ped rvtests.pedigree.ped

############################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">create pedigree for rvtests</entry>
</properties>
"""
}

