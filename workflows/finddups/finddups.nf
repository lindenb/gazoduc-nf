params.paths=""
params.findArgs="-size +100k"
params.max_same = 1000
params.prefix="PREFIX"
params.publishDir="."

workflow {
	ch1 = FIND_FILES(Channel.of(params.paths.split("[\\:]+")).unique())
	ch2 = MERGE(ch1.output.collect())
	ch3 = FINDDUPS(ch2.output.splitText().map{it.trim()} )
	ch4 = MERGE2(ch3.output.collect())
	}

process FIND_FILES {
tag "${dir}"
afterScript "rm -rf TMP"
input:
	val(dir)
output:
	path("by_size.tsv"),emit:output
script:
"""
hostname 1>&2
mkdir -p TMP

find '${dir}' -type f -readable ${params.findArgs}  -printf "%s\t%p\\n" > TMP/jeter1.tsv || echo "cannot find" 1>&2
LC_ALL=C sort -t '\t' -T TMP -k1,1 TMP/jeter1.tsv | uniq > TMP/jeter2.tsv
mv TMP/jeter2.tsv by_size.tsv
"""
}

process MERGE {
tag "N=${L.size()}"
afterScript "rm -rf TMP"
input:
	val(L)
output:
	path("chunks.tsv"),emit:output
script:
"""
mkdir -p TMP OUT
LC_ALL=C sort -t '\t' -T TMP -k1,1 --merge ${L.join(" ")}  | uniq |\
	awk -F '\t' 'BEGIN{P="";N=0;} {if(\$1!=P) {printf("\t%d\\n",N);N=0;} if(N==0) printf("%d",\$1); printf("%s%s",(N==0?"\t":":"),\$2); P=\$1;N++;} END{if(NR>0) printf("\t%d\\n",N);}' |\
	awk -F '\t' '\$3>1 && \$3 < ${params.max_same}' |\
	awk -F '\t' '{L=int(\$1)*1.0;if(L>1E9 || \$3>5) {N=1;} else if(L>1E8) {N=10;} else {N=100;}  printf("%s\t%d\\n",\$0,N);}' > TMP/jeter2.tsv

cut -f 4 TMP/jeter2.tsv | uniq | sort -nr | uniq | while read CSIZE
do
	echo "\${CSIZE}:"
	awk -v S=\${CSIZE} -F '\t' '(\$4==S)' TMP/jeter2.tsv |wc -l
	awk -v S=\${CSIZE} -F '\t' '(\$4==S)' TMP/jeter2.tsv | cut -f1,2 | split -a 9  --lines=\${CSIZE} --additional-suffix=.tsv - "OUT/split.\${CSIZE}."
done

find \${PWD}/OUT/ -type f -name "split.*.tsv" |sort -T TMP > chunks.tsv
test -s chunks.tsv
"""
}

process FINDDUPS {
afterScript "rm -rf TMP"
tag "${file(size2path).name}"
input:
	val(size2path)
output:
	path("output.tsv"),emit:output
script:
"""
mkdir -p TMP
SECONDS=1
touch TMP/output.tsv
cut -f2 "${size2path}" | while IFS=\$"\\t"  read FILES
do
	echo "\${FILES}" | tr ":" "\\n" | sort | uniq | awk '{printf(".\t%s\\n",\$0);}' > TMP/jeter.a
	join -t '\t' -1 1 -2 1 -o '1.2,2.2' TMP/jeter.a TMP/jeter.a | awk -F '\t' '(\$1 < \$2)' | while IFS=\$'\\t' read F1 F2
	do
		echo "\${SECONDS}\${F1}\t\${F2}" >> TMP/log.txt
 		if test -f "\${F1}" && test -f "\${F2}" && cmp --silent "\${F1}" "\${F2}" ; then
			stat  --format '%s\t%U\t%n\t%y\t%Y\t%A' "\${F1}" "\${F2}"  | paste - - >> TMP/output.tsv
		fi
	done
done
LC_ALL=C sort -T TMP -t '\t' -k1,1n TMP/output.tsv > output.tsv
"""
}

process MERGE2 {
tag "N=${L.size()}"
input:
	val(L)
output:
	path("${params.prefix}output.tsv.gz"),emit:output
script:
"""
mkdir -p TMP

echo "size.1 owner.1 file.1 modified.1 timestamp.1 permission.1" > TMP/jeter.a
tr "1" "2" < TMP/jeter.a > TMP/jeter.b

paste TMP/jeter.a TMP/jeter.b | tr " " "\t" > TMP/output.tsv

cat << EOF | tr "\\n" "\\0" | LC_ALL=C sort -T TMP --merge --files0-from=- -t '\t' -k1,1n | uniq >> TMP/output.tsv
${L.join("\n")}
EOF

gzip --best TMP/output.tsv
mv TMP/output.tsv.gz "${params.prefix}output.tsv.gz"
"""
}
