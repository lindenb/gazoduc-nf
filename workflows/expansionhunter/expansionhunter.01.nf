/*

Copyright (c) 2022 Pierre Lindenbaum

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



/** path to indexed fasta reference */
params.reference = ""
params.references = "NO_FILE"
params.prefix = ""
params.publishDir = ""
params.bams = "NO_FILE"
params.help=false
params.extension=1000
params.mapq = 60
params.version="v5.0.0"
params.keepbam=false
params.mode="seeking"
params.catalog="NO_FILE"
params.pedigree="NO_FILE"



include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {parseBoolean;isHg19;getVersionCmd;moduleLoad;runOnComplete} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {SEX_GUESS_02} from '../../subworkflows/sex/sex.guess.02.nf'

workflow {
	EXPANSION_HUNTER_01(params,params.reference,file(params.references),file(params.catalog),file(params.bams),file(params.pedigree))
	}

workflow EXPANSION_HUNTER_01 {
	take:
		meta
		reference
		references
		catalog
		bams
		pedigree
	main:
		version_ch = Channel.empty()
		xhunter = DOWNLOAD_XHUNTER(meta,reference,catalog)
		version_ch = version_ch.mix(xhunter.version)

		ch1 = SEX_GUESS_02(params, reference, references, bams)
		version_ch = version_ch.mix(ch1.version)

		sn_bam = ch1.output.splitCsv(header:true,sep:'\t')

		xh_ch  = RUN_XHUNTER(meta, reference, xhunter.xhunter, xhunter.catalog, sn_bam)
		version_ch = version_ch.mix(xh_ch.version)

		merge_ch = MERGE_VCFS(meta, xh_ch.output.map{T->T[1]}.collect())
		version_ch = version_ch.mix(merge_ch.version)

		version_ch = MERGE_VERSION(meta, "ExpansionHunter", "Expansion hunter", version_ch.collect())

	emit:
		version = version_ch
		pdf  = ch1.pdf
		vcf = merge_ch.vcf
	}

process DOWNLOAD_XHUNTER {
	executor "local"
	afterScript "rm -f jeter.tar.gz"
	input:
		val(meta)
		val(reference)
		path(catalog)
	output:
		path("ExpansionHunter-linux_x86_64/bin/ExpansionHunter"),emit:xhunter
		path("CATALOG/variant_catalog.json"),emit:catalog
		path("version.xml"),emit:version
	script:
		def xversion = meta.expansion_hunter_version?:"5.0.0"
		def url = isHg19(reference) ? "https://github.com/Illumina/RepeatCatalogs/blob/master/hg19/variant_catalog.json?raw=true":""
	"""
	hostname 1>&2
	set -x
	wget -O jeter.tar.gz "https://github.com/Illumina/ExpansionHunter/releases/download/v${xversion}/ExpansionHunter-v${xversion}-linux_x86_64.tar.gz"
	tar xvfz jeter.tar.gz
	rm jeter.tar.gz
	chmod +x "ExpansionHunter-v${xversion}-linux_x86_64/bin/ExpansionHunter"
	# rename directory
	mv "ExpansionHunter-v${xversion}-linux_x86_64"  "ExpansionHunter-linux_x86_64"


	mkdir CATALOG
	if ${!catalog.name.equals("NO_FILE")} ; then
		cp -v "${catalog}" CATALOG/variant_catalog.json
	elif [ ! -z "${url}" ] ; then
		wget -O CATALOG/variant_catalog.json "${url}"
	fi

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">Download expansion hunter and catalog</entry>
	<entry key="expansion.hunter.version">${xversion}</entry>
	<entry key="catalog.url"><url>${url}</url></entry>
	<entry key="version">${getVersionCmd("wget")}</entry>
</properties>
EOF
	"""
	}

process RUN_XHUNTER {
	tag "${row.sample}=${row.sex} / ${file(row.bam).name}"
	afterScript "rm -rf TMP"
	memory "3g"
	cpus 1
	input:
		val(meta)
		val(reference)
		path(executable)
		path(catalog)
		val(row)
	output:
		tuple val(meta),path("${row.new_sample}.vcf.gz"),path("${row.new_sample}.json.gz"),emit:output
		tuple val(meta),path("${row.new_sample}.realigned.cram"),optional:true,emit:bam
		path("version.xml"),emit:version
	script:

	"""
	hostname 1>&2
	${moduleLoad("samtools bcftools")}
        mkdir TMP

	MODE=`grep -o -w -F LocusId '${catalog}' | awk '{N++;} END{printf("%s\\n",(N<1000?"seeking":"streaming"));}'`

	${executable.toRealPath()} \
		--reads "${row.bam}" \
		--reference "${row.reference}" \
		--variant-catalog "${catalog}" \
		--threads ${task.cpus} \
		--output-prefix "TMP/jeter" \
		--region-extension-length ${meta.extension} \
		--analysis-mode "\${MODE}" \
		--sex '${row.sex}' \
		--log-level info 1>&2
	
	if ${!reference.equals(row.reference)} ; then
		TODO
	fi

	# rename sample
	bcftools query -l TMP/jeter.vcf | awk '{printf("%s\t${row.new_sample}\\n",\$1);}' > TMP/jeter.txt
	bcftools reheader --fai "${row.reference}.fai" --samples TMP/jeter.txt TMP/jeter.vcf |\
		awk '/^#CHROM/ {printf("##samples.sex=${row.new_sample}=${row.sex}\\n",S);} {print;}' > TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter.vcf
	rm TMP/jeter.txt
	
	#sort and index
	bcftools annotate --set-id '%INFO/REPID' -O u TMP/jeter.vcf |\
	bcftools sort --max-mem '${task.memory.giga}G' -T TMP -O z -o "TMP/${row.new_sample}.vcf.gz
	bcftools index -t "TMP/${row.new_sample}.vcf.gz"

	# save json
	gzip --best TMP/jeter.json 

	if  ${parseBoolean(meta.keepbam)} ; then
		samtools sort -T TMP/sort  -@ ${task.cpus} -O CRAM --reference "${row.reference}" -o "${row.new_sample}.realigned.cram" "TMP/jeter_realigned.bam"
		samtools index "${row.new_sample}.realigned.cram"
	fi


	mv TMP/${row.new_sample}.vcf.gz ./
	mv TMP/${row.new_sample}.vcf.gz.tbi ./
	mv TMP/jeter.json.gz "${row.new_sample}.json.gz"

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">run expansion hunter</entry>
	<entry key="mode">\${MODE}</entry>
	<entry key="sample">${row.new_sample}</entry>
	<entry key="bam">${row.bam}</entry>
	<entry key="version">${getVersionCmd("samtools bcftools")}</entry>
</properties>
EOF
	"""
	}


process MERGE_VCFS {
	tag "N=${L.size()}"
	cpus 10
	input:
		val(meta)
		val(L)
	output:
		path("${meta.prefix?:""}merged.bcf"),emit:vcf
		path("vcfs.list"),emit:vcfs_list
		path("version.xml"),emit:version
	script:
	"""
	hostname 1>&2
	${moduleLoad("bcftools")}
	set -o pipefail

cat << EOF > vcfs.list
${L.join("\n")}
EOF

	bcftools merge --threads ${task.cpus} -m id -O u --file-list vcfs.list |\
		bcftools sort  -T . -o "${meta.prefix?:""}merged.bcf" -O b

	bcftools index "${meta.prefix?:""}merged.bcf"

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">merge vcfs</entry>
	<entry key="vcf.count">${L.size()}</entry>
	<entry key="version">${getVersionCmd("bcftools")}</entry>
</properties>
EOF
	"""
	}


process miniKit {
	executor "local"
	output:
		path("minikit.jar") into minikitjar
	when:
		true
	script:
	"""
cat << __EOF__ > MiniKit.java
import java.io.*;
import java.nio.file.*;
import java.nio.charset.*;
import java.util.*;
import java.util.stream.*;
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.*;
import htsjdk.variant.vcf.*;
import org.broadinstitute.hellbender.utils.FisherExactTest;

public class MiniKit {
private final Set<String> cases = new HashSet<>();
private final Set<String> controls = new HashSet<>();

private int allele2length(final Allele a) {
	if(!a.isSymbolic()) return -1;
	String s = a.getDisplayString();
	if(!s.startsWith("<STR")) return -1;
	if(!s.endsWith(">")) return -1;
	return Integer.parseInt(s.substring(4,s.length()-1));
	}

void instanceMainWithExit(final String args[]) {
	try {
		try(BufferedReader br = Files.newBufferedReader(Paths.get(args[0]),Charset.defaultCharset()) ) {
			br.lines().map(L->L.split("[\\t]")).forEach(T->{
				if(T[1].equals("case")) cases.add(T[0]);
				if(T[1].equals("control")) controls.add(T[0]);
				});
			}
		catch(final Throwable err) {
			err.printStackTrace();
			System.exit(-1);
			}
		if(cases.isEmpty() || controls.isEmpty()) {
			System.err.println("input error");
			System.exit(-1);
			}
		System.err.println(""+cases.size()+" "+controls.size());
		try(final VCFIterator r = new VCFIteratorBuilder().open(System.in)) {
		final VCFHeader h = r.getHeader();
		final int  array[][] = new int[2][];
		array[0] = new int[2];
		array[1] = new int[2];
		while(r.hasNext()) {
			final VariantContext ctx = r.next();
			final int ref_num = ctx.getAttributeAsInt("REF",-1);
			final Set<Integer> lengths = new HashSet<>(ctx.getAlternateAlleles().
					stream().
					map(A->allele2length(A)).
					filter(L->L.intValue() > ref_num ).
					collect(Collectors.toList())
					);
			if(lengths.isEmpty()) continue;
			double best = 1.0;
			String line="";
			for(Integer len : lengths) {
				final int curr_len = len;
				int case_alt = 0;
				int case_ref = 0;
				int ctrl_alt = 0;
				int ctrl_ref = 0;
				for(final Genotype gt : ctx.getGenotypes()) {
					if(gt.isNoCall()) continue;
					final String sn = gt.getSampleName();
					if(!(controls.contains(sn) || cases.contains(sn))) continue;
					final boolean is_control = controls.contains(sn);
					final boolean has_len = gt.getAlleles().stream().mapToInt(A->allele2length(A)).anyMatch(L-> L >= curr_len);

					     if( is_control &&  has_len) { ctrl_alt++;}
					else if( is_control && !has_len) { ctrl_ref++;}
					else if(!is_control &&  has_len) { case_alt++;}
					else if(!is_control && !has_len) { case_ref++;}
					}
				array[0][0] = ctrl_alt;
				array[0][1] = ctrl_ref;
				array[1][0] = case_alt;
				array[1][1] = case_ref;
				final double fisher= FisherExactTest.twoSidedPValue(array);
				if(fisher > best )  continue;
				best = fisher;
				line = 	 ctx.getAttributeAsString("VARID",".") + "\t" +
						ctx.getContig() + "\t" + ctx.getStart() +
						"\t" + fisher +
						"\t<STR" + (curr_len) + ">\t" +
						ctx.getAttributeAsString("RU",".") + "\t" +
						ctrl_alt + "\t" + ctrl_ref + "\t" + case_alt + "\t" + case_ref
						;
				}
			if(line.isEmpty()) continue;
			System.out.println(line);
			}
		}
		}
	catch(final Throwable err) {
		err.printStackTrace();
		System.exit(-1);
		}
	}

public static void main(final String args[]) {
	new MiniKit().instanceMainWithExit(args);
	}

}
__EOF__

mkdir tmp
javac -cp ${params.classpath} -d tmp MiniKit.java
jar cvf minikit.jar -C tmp .

"""
}

process fisherTest {
	input:
		val(sex) from all_sex
		val(jar) from minikitjar
		val(vcf) from merged_bcf
	output:
		path("${params.prefix}fisher.tsv.gz") into fisher_tsv
	script:
	"""
	module load bcftools/0.0.0
	set -o pipefail
	bcftools view "${vcf}" |\
		java -cp ${params.classpath}:${jar} MiniKit "${sex}" > jeter.tsv

	LC_ALL=C sort -T . -t '\t' -k4,4g jeter.tsv > "${params.prefix}fisher.tsv"
	gzip --best "${params.prefix}fisher.tsv"
	rm jeter.tsv
	"""
}

process installPackages {
output:
	path("LIB") into r_packages_lib
script:
"""
module load r/3.6.3
hostname 1>&2
mkdir LIB
cat << EOF | R --vanilla
install.packages(c("qqman"),repos="http://cran.r-project.org",lib="\${PWD}/LIB")
EOF
"""
}


process plotIt {
afterScript "rm -f jeter.tsv jeter.R"
input:
	path rlib from r_packages_lib
	path(fisher) from fisher_tsv
output:
	path("paths.txt") into assoc_plot	
script:
"""
module load r/3.6.3
hostname 1>&2

echo "SNP\tCHR\tBP\tP" > jeter.tsv
gunzip -c "${fisher}" | cut -f 1,2,3,4 | awk -F '\t' '{OFS="\t";C=\$2;gsub(/^chr/,"",C);if(C=="X") C="23"; if(C=="Y") C="24";\$2=C;print;}' >> jeter.tsv


cat << '__EOF__' > jeter.R
library("qqman",lib.loc="${rlib.toRealPath()}")
T1 <- read.table("jeter.tsv",header=TRUE,sep="\t",stringsAsFactors=FALSE)

if(nrow(T1)>0) {
png("${params.prefix}fisher.manhattan.png")
manhattan(T1,main="${params.prefix}fisher");
dev.off()

png("${params.prefix}fisher.qqplot.png")
qq(T1\$P,main="${params.prefix}fisher");
dev.off()
}
__EOF__

R --vanilla < jeter.R || true

find \${PWD} -type f -name "*.png" >> paths.txt
"""
}



process zipBcf {
publishDir "${params.publishDir}" ,  mode: 'copy', overwrite: true
input:
	val(assoc) from assoc_plot
	val(fisher) from fisher_tsv
	val(readme) from readme_md
	val(sex) from all_sex
	val(merged) from merged_bcf
	val(vcfspath) from vcfs_list
output:
        path("${params.prefix}all.zip") into merged
when:
     true
script:
"""

#cp "${vcfspath}" jeter.txt

echo "${sex}" >> jeter.txt
echo "${merged}" >> jeter.txt
echo "${assoc}" >> jeter.txt
echo "${fisher}" >> jeter.txt


zip -j -@ -0 "${params.prefix}all.zip" < jeter.txt

rm jeter.txt
"""
}


	
runOnComplete(workflow);

