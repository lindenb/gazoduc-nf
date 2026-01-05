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
nextflow.enable.dsl=2

include {isHg19;runOnComplete;moduleLoad;getKeyValue;hasFeature;getVersionCmd;parseBoolean} from '../../../modules/utils/functions.nf'
include {BURDEN_SAMPLE_WGSELECT_PART_01}  from '../../../subworkflows/burden/burden.samples.wgselect.part.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.nf'
include {RVTESTS_REHEADER_01} from '../../../modules/rvtests/rvtests.reheader.01.nf'
include {RVTESTS_POST_PROCESS} from '../../../subworkflows/rvtests/rvtests.post.process.01.nf'
include {CONCAT_FILES_01} from '../../../modules/utils/concat.files.nf'

params.reference=""
params.pedigree=""
params.vcf=""
params.disableFeatures="";
params.help=false
/** use only introns overlapping this bed */
params.bed= "NO_FILE"
/** in  intron, only keep/annotate in this interval(s) */
params.bed_cluster_method = " --size 1mb "
params.soacn = "" /* empty, take all */
params.with_utr5=true
params.with_utr3=true

if(params.help) {
  log.info"""
## About

Burden for 1st intron.

## Author

Pierre Lindenbaum

## Options

  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --vcf <file> path to a indexed VCF or BCF file. If file ends with '.list' is a list of path to one VCF per contig [REQUIRED]
  * --pedigree <file> jvarkit formatted pedigree. phenotype MUST be case|control. Sex MUST be male|female|unknown
  * --bed <file> optional bed file to limit the analysis to the genes overlapping a  bed file.
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume workflow.nf \\
        --publishDir output \\
        --prefix "analysis." \\
        --reference /path/to/reference.fasta \\
        --vcf /path/to/my.vcf.gz \\
        --pedigree /path/to/input.ped \
```

## Workflow

![workflow](./workflow.svg)

"""
exit 0
}

workflow {
		burden_ch = BURDEN_UTR(params, params.reference, params.vcf, file(params.pedigree),file(params.bed))
		ZIPIT(params,burden_ch.zip.collect())
		}

workflow BURDEN_UTR {
	take:
		meta
		reference
		vcf
		pedigree
		bed
	main:

		version_ch = Channel.empty()
		to_zip = Channel.empty()
		

		utr_ch = EXTRACT_UTR(meta, reference, bed)
		version_ch = version_ch.mix(utr_ch.version)
	
		ch1_ch = BURDEN_SAMPLE_WGSELECT_PART_01(meta,reference,vcf, pedigree, utr_ch.merged_bed)
		version_ch = version_ch.mix(ch1_ch.version)

		each_setfilelist_ch = utr_ch.output.splitText().
				map{S->file(S.trim())}
		
		header_ch = RVTESTS_REHEADER_01(meta, reference)
		version_ch = version_ch.mix(header_ch.version)

		assoc_ch = RVTEST_UTR(meta, reference, ch1_ch.contig_vcfs, ch1_ch.rvtest_pedigree, header_ch.output, each_setfilelist_ch)
		version_ch = version_ch.mix(assoc_ch.version)
		
		concat_ch = CONCAT_FILES_01(meta,assoc_ch.output.collect())
		version_ch = version_ch.mix(concat_ch.version)

		digest_ch = RVTESTS_POST_PROCESS(meta, reference, file(vcf).name ,concat_ch.output)
                version_ch = version_ch.mix(digest_ch.version)
		to_zip = to_zip.mix(digest_ch.zip)

		version_ch = MERGE_VERSION(meta, "burden UTR", "Burden UTR ${vcf}", version_ch.collect())
		to_zip = to_zip.mix(version_ch.version)

		html = VERSION_TO_HTML(params,version_ch.version)
		to_zip = to_zip.mix(html.html)
		
	emit:
		version = version_ch
		zip = to_zip
	}


process EXTRACT_UTR {
memory "3g"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
	path(bed)
output:
	path("merged.bed"),emit:merged_bed
	path("setfile.list"),emit:output
	path("version.xml"),emit:version
script:
	def url = isHg19(reference)?"http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/ensGene.txt.gz":""
"""
hostname 1>&2
${moduleLoad("java bedtools htslib")}
set -o pipefail

test ! -z "${url}"

mkdir -p TMP

cat << EOF > TMP/Minikit.java
import java.util.*;
import java.io.*;
public class Minikit {
	private static class Interval {
		final int start;
		final int end;
		Interval(int start,int end) {
			this.start = start;
			this.end = end;
			}
		}
	public static void main(String[] args) {
	final boolean with_utr5 = ${parseBoolean(meta.with_utr5)};
	final boolean with_utr3 = ${parseBoolean(meta.with_utr3)};
	try(BufferedReader br=new BufferedReader(new InputStreamReader(System.in))) {
		String line;
		while((line=br.readLine())!=null) {
			final String[] tokens = line.split("[\t]");
			final int cdsStart = Integer.parseInt(tokens[6]);
			final int cdsEnd = Integer.parseInt(tokens[7]);
			if(cdsStart>= cdsEnd) continue;

			final int nExons = Integer.parseInt(tokens[8]);
			final int[] exonStarts = Arrays.stream(tokens[ 9].split("[,]")).mapToInt(S->Integer.parseInt(S)).toArray();
			final int[] exonEnds   = Arrays.stream(tokens[10].split("[,]")).mapToInt(S->Integer.parseInt(S)).toArray();
			final char strand = tokens[3].charAt(0);
			final List<Interval> L5 = new ArrayList<>();
			final List<Interval> L3 = new ArrayList<>();
			if((strand=='+' && with_utr5) || (strand=='-' && with_utr3)) {
				final List<Interval> L = (strand=='+'?L5:L3);
				for(int i=0;i<nExons;i++) {
					if(exonStarts[i] >= cdsStart) break;
					L.add(new Interval(exonStarts[i],Math.min(cdsStart,exonEnds[i])));
					}
				}
			if((strand=='-' && with_utr5) || (strand=='+' && with_utr3)) {
				final List<Interval> L = (strand=='+'?L3:L5);
				for(int i=nExons-1;i>=0;i--) {
					if(exonEnds[i] <= cdsEnd) break;
					L.add(new Interval(Math.max(cdsEnd,exonStarts[i]),exonEnds[i]));
					}
				}
			final List<Interval> L53 = new ArrayList<>(L5);
			L53.addAll(L3);

			for(int side=0;side<3;++side) {
				final List<Interval> L;
				final String suffix;
				switch(side) {
					case 0: L=L5; suffix="UTR_5"; break;
					case 1: L=L3; suffix="UTR_3"; break;
					default: L=L53; suffix="UTR_ALL"; break;
					}
				if(L.isEmpty()) continue;
				Collections.sort(L,(A,B)->Integer.compare(A.start,B.start));
				System.out.print(tokens[12]+"_"+tokens[1]+"_"+suffix+"\t");
				for(int i=0;i< L.size();i++) {
					if(i>0)  System.out.print(",");
					System.out.print(tokens[2]);
					System.out.print(":");
					System.out.print(L.get(i).start+1);
					System.out.print("-");
					System.out.print(L.get(i).end);
					}
				System.out.println();
				}
			}
		}
	catch(Throwable err) {
		err.printStackTrace();
		System.exit(-1);
		}
	}
}
EOF

javac -d TMP -sourcepath TMP TMP/Minikit.java

wget -O - "${url}" | gunzip -c |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 3 --convert SKIP |\
	java -cp TMP Minikit |\
	sort -T TMP -t '\t' -k2,2 --unique > TMP/jeter.setfile


if ${!bed.name.equals("NO_FILE")} ; then
	${bed.name.endsWith(".gz")?"gunzip -c ":"cat"} "${bed}" |\
		cut -f1,2,3 |
		java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP |\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\
		bedtools merge  > TMP/jeter.bed
	bgzip TMP/jeter.bed
	tabix -p bed --force TMP/jeter.bed.gz

	java -jar \${JVARKIT_DIST}/setfiletools.jar -R "${reference}" intersectbed TMP/jeter.bed.gz TMP/jeter.setfile > TMP/jeter2.setfile
	mv TMP/jeter2.setfile TMP/jeter.setfile

fi

# group the set files by pool of 'x' bases
mkdir BEDS
java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/setfiletools.jar -R "${reference}" cluster --out BEDS  --size "1Mb" TMP/jeter.setfile

find \${PWD}/BEDS -type f -name "*.setfile" > setfile.list
test -s setfile.list

# create merged bed

java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/setfiletools.jar -R "${reference}" tobed TMP/jeter.setfile |\
	cut -f 1,2,3 |\
	sort -T TMP -t '\t' -k1,1 -k2,2n |\
	bedtools merge > TMP/jeter.bed

mv TMP/jeter.bed  merged.bed

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">extract UTR from UCSC</entry>
        <entry key="utr5">${meta.with_utr5}</entry>
        <entry key="utr3">${meta.with_utr3}</entry>
        <entry key="description">extract UTR from UCSC</entry>
	<entry key="url"><url>${url}</url></entry>
	<entry key="versions">${getVersionCmd("bedtools tabix jvarkit/bedrenamechr jvarkit/setfiletools")}</entry>
	<entry key="bed">${bed}</entry>
</properties>
EOF
"""
}


process RVTEST_UTR {
tag "${setfile.name}"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
	path(vcfs)
	path(pedigree)
	path(reheader)
	path(setfile)
output:
	path("assoc.list"),emit:output
	path("version.xml"),emit:version
script:
	def  rvtest_params = "--burden 'cmc,exactCMC,zeggini' --kernel 'skato'"

	// --burden cmc,zeggini,mb,fp,exactCMC,cmcWald,rarecover,cmat --vt price,analytic --kernel 'skat[nPerm=1000],kbac,skato'

"""
hostname 1>&2
${moduleLoad("rvtests jvarkit bedtools bcftools")}
set -o pipefail
set -x
mkdir -p TMP ASSOC


java -jar \${JVARKIT_DIST}/setfiletools.jar -R "${reference}" tobed "${setfile}" |\
	cut -f 1,2,3 | sort -T TMP -t '\t' -k1,1 -k2,2n | bedtools merge > TMP/jeter.bed

#remove chr prefix
java -jar \${JVARKIT_DIST}/setfiletools.jar -R "${reference}" view --trim-chr '${setfile}' > TMP/jeter.setfile

bcftools concat -a --regions-file TMP/jeter.bed --file-list "${vcfs}" -O u |\
	bcftools annotate  --rename-chrs "${reheader}" -O z -o TMP/jeter.vcf.gz
bcftools index -t TMP/jeter.vcf.gz


rvtest  --noweb \
        --inVcf TMP/jeter.vcf.gz \
	--setFile TMP/jeter.setfile \
	--pheno "${pedigree}" \
	--out "ASSOC/part" \
	${rvtest_params} 1>&2 2> TMP/last.rvtest.log


find \${PWD}/ASSOC -type f -name "part*assoc" > assoc.list

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">invoke rvtest for setfile</entry>
	<entry key="rvtest.path">\$(which rvtest)</entry>
	<entry key="versions">${getVersionCmd("bedtools rvtest bcftools jvarkit/setfiletools")}</entry>
	<entry key="setfile">${setfile}</entry>
</properties>
EOF
"""
}


process ZIPIT {
tag "N=${files.size()}"
publishDir "${meta.publishDir}" , mode: 'copy', overwrite: true
input:
	val(meta)
        val(files)
output:
        path("${meta.prefix?:""}archive.zip")
when:
        !meta.getOrDefault("publishDir","").trim().isEmpty()
script:
        prefix = meta.getOrDefault("prefix","")
"""

mkdir "${prefix}archive"

cat << EOF | while read F ; do ln -s "\${F}" "./${prefix}archive/" ; done
${files.join("\n")}
EOF

zip -r "${prefix}archive.zip" "${prefix}archive"
"""
}


runOnComplete(workflow);
