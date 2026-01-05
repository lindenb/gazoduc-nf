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
include {GS_SIMPLE_01 as GS1; GS_SIMPLE_01 as GS2} from '../../../modules/gs/gs.simple.01.nf'


params.reference=""
params.pedigree=""
params.vcf=""
params.help=false
params.disableFeatures=""
params.keyRegex=""

if(params.help) {
  log.info"""
## About

optimize burden using sliding window of variants.

## Author

Pierre Lindenbaum

## Options

  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --vcf <file> path to a indexed VCF or BCF file. If file ends with '.list' is a list of path to one VCF per contig [REQUIRED]
  * --pedigree <file> jvarkit formatted pedigree. phenotype MUST be case|control. Sex MUST be male|female|unknown
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
		burden_ch = OPTIMIZE_BURDEN(params, params.reference, params.vcf, file(params.pedigree))
		//ZIPIT(params,burden_ch.zip.collect())
		}

workflow OPTIMIZE_BURDEN {
	take:
		meta
		reference
		vcf
		pedigree
	main:

		version_ch = Channel.empty()
		to_zip = Channel.empty()
		
		vcfbed_ch = VCF2BED(meta,reference,vcf)
		version_ch = version_ch.mix(vcfbed_ch.version)
	
		ch1_ch = BURDEN_SAMPLE_WGSELECT_PART_01(meta,reference,vcf, pedigree, vcfbed_ch.bed)
		version_ch = version_ch.mix(ch1_ch.version)

		
		header_ch = RVTESTS_REHEADER_01(meta, reference)
		version_ch = version_ch.mix(header_ch.version)

		genes_ch = SPLIT_GENES(meta,reference,header_ch.output,ch1_ch.contig_vcfs)
		version_ch = version_ch.mix(genes_ch.version)

		rows = genes_ch.output.splitCsv(header:true,sep:'\t').
			filter{T->(meta.keyRegex.isEmpty() || T.key.matches(meta.keyRegex))}.
			filter{T->(T.Count_Variants as int)>2}

		window_ch = WINDOW_GENE(meta, reference, ch1_ch.rvtest_pedigree, rows )
		version_ch = version_ch.mix(window_ch.version)

		assoc_ch = RVTEST_SETFILE(meta, reference, window_ch.output.splitCsv(header:true,sep:','))
		version_ch = version_ch.mix(assoc_ch.version)

		ch3_ch = assoc_ch.output.flatMap{T->[ [[T[0],"CMCFisherExact"],T[1]] , [[T[0],"CMC"],T[2] ] ]}
		plot_ch = PLOT_IT(meta, ch3_ch.groupTuple())
		version_ch = version_ch.mix(plot_ch.version)

		pdf0_ch = GS1(meta.plus("with_header":true),plot_ch.pdf.groupTuple())
		version_ch = version_ch.mix(pdf0_ch.version)

		pdf_ch = GS2(meta,pdf0_ch.output.map{T->["optimize",T ]}.groupTuple())
		version_ch = version_ch.mix(pdf_ch.version)


/*

		version_ch = MERGE_VERSION(meta, "burden UTR", "Burden UTR ${vcf}", version_ch.collect())
		to_zip = to_zip.mix(version_ch.version)

		html = VERSION_TO_HTML(params,version_ch.version)
		to_zip = to_zip.mix(html.html)
		
	emit:
		version = version_ch
		zip = to_zip
*/
	}


process VCF2BED {
tag "${file(vcf).name}"
executor "local"
input:
	val(meta)
	val(reference)
	val(vcf)
output:
	path("vcf.bed"),emit:bed
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bcftools bedtools")}
set -o pipefail

bcftools query -f '%CHROM\t%POS0\t%END\\n' "${vcf}" |\
	sort -T . -t '\t' -k1,1 -k2,2n |\
	bedtools merge > vcf.bed

# check one contig
test 1 -eq `cut -f 1 vcf.bed | sort | uniq | wc -l`

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">BED FROM VCF</entry>
	<entry key="versions">${getVersionCmd("bedtools bcftools")}</entry>
	<entry key="vcf">${vcf}</entry>
</properties>
EOF
"""
}


process SPLIT_GENES {
tag "${vcfs.name}"
memory "3g"
input:
	val(meta)
	val(reference)
	path(reheader)
	path(vcfs)
output:
	path("splitgene.mf"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("jvarkit bcftools")}
set -o pipefail
set -x

mkdir -p TMP
bcftools concat --file-list "${vcfs}" -a --remove-duplicates -O u |\
	bcftools annotate -O v --rename-chrs "${reheader}" |\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/vcfgenesplitter.jar \
		-E 'ANN/GeneId ANN/FeatureId VEP/GeneId VEP/Ensp VEP/Feature' \
		--manifest splitgene.mf \
		-o \${PWD}/TMP
		
find  TMP  -type f -name "*.vcf.gz" | while read F
do
	bcftools index --force --tbi "\${F}"
done

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">SPLIT VCF PER GENE</entry>
	<entry key="versions">${getVersionCmd("jvarkit/vcfgenesplitter bcftools")}</entry>
	<entry key="vcfs">${vcfs}</entry>
</properties>
EOF
"""
}

process WINDOW_GENE {
tag "${file(row.path).name}"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
	path(pedigree)
	val(row)
output:
	path("output.csv"),emit:output
	path("version.xml"),emit:version
script:
	def name= row.splitter.replaceAll("/","_")+"_"+row.gene+"_"+row.key;
	def vcf = row.path
"""
hostname 1>&2
${moduleLoad("bcftools picard")}
set -o pipefail
set -x
mkdir -p TMP BEDS

cat << EOF > TMP/Minikit.java
import java.util.*;
import java.io.*;
import java.util.function.Function;
import java.util.stream.Collectors;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;
import htsjdk.samtools.util.*;

public class Minikit {
private int instanceMain(final String[] args) {
	try {
		final List<Interval> intervals = new ArrayList<>();

		try(VCFIterator iter=new VCFIteratorBuilder().open(System.in)) {
			while(iter.hasNext()) {
				final VariantContext ctx = iter.next();
				intervals.add(new Interval(ctx.getContig(),ctx.getStart(),ctx.getEnd()));
				}
			}
 		Collections.sort(intervals);
		int i=0;
		while(i+1<intervals.size()) {
			if(intervals.get(i).intersects(intervals.get(i+1))) {
				intervals.set(i,new Interval(
					intervals.get(i).getContig(),
					Math.min(intervals.get(i).getStart(), intervals.get(i+1).getStart()),
					Math.max(intervals.get(i).getEnd(), intervals.get(i+1).getEnd())
					));
				intervals.remove(i+1);
				}
			else
				{
				i++;
				}
			}
		System.err.println("N="+intervals.size());
		for(int n=2; n<= intervals.size();n++) {
			try(PrintWriter pw = new PrintWriter("BEDS/"+n+".setfile")) {
			for(int x=0;(x+n) <= intervals.size(); ++x) {
				pw.println("${name}_b"+x+"_n"+n+"\t"+
					intervals.subList(x,x+n).
					stream().
					map(R->R.getContig()+":"+R.getStart()+"-"+R.getEnd()).
					collect(Collectors.joining(","))
					);
				}
			pw.flush();
			}
			}
		return 0;
		}
	catch(Throwable err) {
		err.printStackTrace();
		return -1;
		}
	}

private void instanceMainWithExit(final String[] args) {
	System.exit(instanceMain(args));
	}

public static void main(String[] args) {
	new Minikit().instanceMainWithExit(args);
	}
}
EOF

javac -cp "\${PICARD_JAR}" -d TMP -sourcepath TMP TMP/Minikit.java

bcftools view -G "${vcf}" | java -cp "\${PICARD_JAR}:TMP" Minikit > TMP/input.setfile

find \${PWD}/BEDS -type f -name "*.setfile" |\
awk 'BEGIN{printf("vcf,pedigree,name,setfile\\n");} {printf("${vcf},${pedigree.toRealPath()},${name},%s\\n",\$0);}' > output.csv

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">split into setfile</entry>
	<entry key="versions">${getVersionCmd("awk javac")}</entry>
	<entry key="name">${name}</entry>
	<entry key="vcf">${vcf}</entry>
</properties>
EOF
"""	
}





process RVTEST_SETFILE {
tag "${row.name} ${file(row.setfile).name}"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
	val(row)
output:
	tuple val("${row.name}"),path("ASSOC/part.CMCFisherExact.assoc"),path("ASSOC/part.CMC.assoc"),emit:output
	path("version.xml"),emit:version
script:
	def name= row.name
	def vcf = row.vcf
	def  rvtest_params = "--burden 'cmc,exactCMC' "
"""
hostname 1>&2
${moduleLoad("rvtests")}
set -o pipefail
set -x
mkdir -p ASSOC TMP


rvtest  --noweb \
        --inVcf "${vcf}" \
	--setFile ${row.setfile} \
	--pheno "${row.pedigree}" \
	--out "ASSOC/part" \
	${rvtest_params} 1>&2 2> TMP/last.rvtest.log

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">invoke rvtest for setfile</entry>
	<entry key="rvtest.path">\$(which rvtest)</entry>
	<entry key="versions">${getVersionCmd("rvtest")}</entry>
	<entry key="name">${name}</entry>
	<entry key="vcf">${vcf}</entry>
	<entry key="setfile">${row.setfile}</entry>
</properties>
EOF
"""	
}


process PLOT_IT {
tag "${key[0]} ${key[1]} N=${L.size()}"
afterScript "rm -rf TMP"
input:
	val(meta)
	tuple val(key),val(L)
output:
	tuple val("${key[1]}"),path("${meta.prefix?:""}${key[1]}.${key[0]}.pdf"),emit:pdf
	path("version.xml"),emit:version
script:
	def name  = key[0];
	def type = key[1];
	def column = type.equals("CMC")?"7":"10"
"""
hostname 1>&2
${moduleLoad("r")}

mkdir -p TMP
export LC_ALL=C

cat ${L.join(" ")} |\
	awk '/^Range/ {next;} {N=split(\$2,a,/[,:-]/);for(i=1;i<=N;i+=3) {B=int(a[i+1]);E=int(a[i+2]);if(i==1 || B<m) m=B; if(i==1 || E>M) M=E;} printf("%s\t%s\t%s\\n",m,M,\$${column});}' | sort -T . | uniq > jeter.tsv

cat << "__EOF__" > TMP/jeter.R
T1<-read.table("jeter.tsv",sep="\\t",header=FALSE,col.names=c("start","end","pvalue"),colClasses=c("integer","integer","numeric"))
head(T1)
pdf("${meta.prefix?:""}${key[1]}.${key[0]}.pdf")
color <- rgb(0.4, 0.4, 1.0, 0.2)
plot(1, type="n", main="${type} / ${name} /",xlab="position", ylab="-log10(pvalue)", xlim=c(min(T1\$start),max(T1\$end)), ylim=c(0, -1.0 * log10(min(T1\$pvalue))))
for(i in 1:nrow(T1)) {
	y <- log10(T1[i,]\$pvalue) * -1.0
	segments(x0=T1[i,]\$start,x1=T1[i,]\$end,y0=y,y1=y,col=color)
}
dev.off()
__EOF__

R --vanilla < TMP/jeter.R

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">plot rvtest</entry>
        <entry key="type">${key[1]}</entry>
	<entry key="versions">${getVersionCmd("R awk javac")}</entry>
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
