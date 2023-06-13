
workflow VG_CONSTRUCT_01 {
take:
	meta
	reference
	vg_exec
	vcf
	bed
main:
	version_ch = Channel.empty()
	if(bed.name.equals("NO_FILE")) {
		each_bed = Channel.fromPath(reference+".fai").splitCsv(header:false,sep:'\t').map{T->[T[0],1,T[1]]}
		}
	else
		{
		each_bed = Channel.fromPath(bed).splitCsv(header:false,sep:'\t').map{T->[T[0],T[1],T[2]]}.take(2)
		}
	construct_ch = CONSTRUCT([:],reference, vg_exec, vcf, each_bed)
	merge_ch = MERGE([:], reference, vg_exec, construct_ch.output.collect())

	vg = merge_ch.output.map{T->[vg:T[0],gcsa:T[1],xg:T[2],snarls:T[3],reference:reference]}
emit:
	version=version_ch
	vg = vg
}

process CONSTRUCT {
tag "${contig}:${start}-${end}"
cpus 1
memory "3g"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
	val(vg_exec)
	val(vcf)
	tuple val(contig),val(start),val(end)
output:
	path("path.txt"),emit:output
script:
"""
hostname 1>&2
mkdir -p TMP
export TMPDIR="\${PWD}/TMP"

${vg_exec} construct \
	--reference "${reference}" \
	--vcf "${vcf}" \
	--region "${contig}:${start}-${end}" \
	--threads "${task.cpus}" \
	--handle-sv > TMP/jeter.vg

mv -v TMP/jeter.vg "${contig}_${start}_${end}.vg"
echo "\${PWD}/${contig}_${start}_${end}.vg" > path.txt
"""
}

process MERGE {
tag "N=${L.size()}"
cpus 1
memory "10g"
input:
	val(meta)
	val(reference)
	val(vg_exec)
	val(L)
output:
	tuple path("${file(reference).getSimpleName()}.vg"),
	path("${file(reference).getSimpleName()}.gcsa"),
	path("${file(reference).getSimpleName()}.xg"),
	path("${file(reference).getSimpleName()}.snarls"),emit:output
script:
	def name = file(reference).getSimpleName()
"""
hostname 1>&2
mkdir -p TMP
export TMPDIR="\${PWD}/TMP"

cat ${L.join(" ")} > TMP/jeter.list

## https://github.com/vgteam/vg/issues/2515
${vg_exec} ids -j `cat TMP/jeter.list`

xargs -a TMP/jeter.list cat > TMP/jeter.vg.tmp
mv TMP/jeter.vg.tmp TMP/jeter.vg



${vg_exec} prune -M32 --restore-paths  TMP/jeter.vg > TMP/pruned.vg
${vg_exec} index --progress --temp-dir TMP --gcsa-out TMP/jeter.gcsa TMP/pruned.vg 
${vg_exec} index --progress --temp-dir TMP --xg-name  TMP/jeter.xg TMP/pruned.vg
${vg_exec} snarls TMP/pruned.vg > TMP/ref.snarls

mv TMP/ref.snarls "${name}.snarls"
mv TMP/pruned.vg "${name}.vg"
mv TMP/jeter.gcsa "${name}.gcsa"
mv TMP/jeter.gcsa.lcp "${name}.gcsa.lcp"
mv TMP/jeter.xg "${name}.xg"
"""
}
