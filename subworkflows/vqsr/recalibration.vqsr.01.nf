/*

Copyright (c) 2023 Pierre Lindenbaum

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
include {moduleLoad;escapeXml} from '../../modules/utils/functions.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include { PICARD_GATHER_VCFS_01 } from '../../modules/picard/picard.gathervcfs.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'





workflow VARIANT_VQSR_01 {
	take:
		meta
		genomeId
		vcf
	main:
		version_ch = Channel.empty()


		vcf2bed_ch = BCF_TO_VCF([:],vcf)
		version_ch = version_ch.mix(vcf2bed_ch.version)

		compile_ch= COMPILE_MINIKIT([:])
		version_ch = version_ch.mix(compile_ch.version)

		contigs_ch = vcf2bed_ch.contigs.splitCsv(header:false,sep:'\t').map{T->T[0]}

		recal_snp_ch = RECALIBRATE_SNP([:], genomeId, vcf2bed_ch.vcf,compile_ch.jar)
		version_ch = version_ch.mix(recal_snp_ch.version)

		snp2_ch = recal_snp_ch.output.map{T->[
			"vcf":T[0],
			"recal_snp_vcf":T[1],
			"recal_snp_tranches":T[2]
			]}

		
		recal_indel_ch = RECALIBRATE_INDEL([:],genomeId, vcf2bed_ch.vcf,compile_ch.jar)
		version_ch = version_ch.mix(recal_indel_ch.version)



		indel2_ch = recal_indel_ch.output.map{T->[
			"vcf":T[0],
			"recal_indel_vcf":T[1],
			"recal_indel_tranches":T[2]
			]}

		ch1_ch = contigs_ch.combine(snp2_ch).
			map{T->T[1].plus("contig":T[0])}.
			combine(indel2_ch).
			map{T->T[0].plus(T[1])}

		apply_snp_ch = APPLY_RECALIBRATION_SNP([:], genomeId ,ch1_ch)
		version_ch = version_ch.mix(apply_snp_ch.version)

		// the vcf is now the vcf generated by apply_recal_snp
		ch2_ch = apply_snp_ch.output.map{T->T[0].plus("vcf":T[1])}

		apply_indel_ch = APPLY_RECALIBRATION_INDEL([:], genomeId ,ch2_ch)
		version_ch = version_ch.mix(apply_indel_ch.version)

		cat_files_ch = COLLECT_TO_FILE_01([:], apply_indel_ch.vcf.collect())
		version_ch = version_ch.mix(cat_files_ch.version)

		gather_ch = PICARD_GATHER_VCFS_01(["suffix":".bcf"], cat_files_ch.output)
		version_ch = version_ch.mix(gather_ch.version)
		
		version_ch = MERGE_VERSION("Variant Recalibration", version_ch.collect())

	emit:
		version = version_ch
		vcf = gather_ch.vcf
		index = gather_ch.index
	}


process BCF_TO_VCF {
input:
	val(meta)
	val(vcf)
output:
	path("reformat.vcf.gz"),emit:vcf
	path("reformat.vcf.gz.tbi"),emit:index
	path("chroms.txt"),emit:contigs
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}
set -o pipefail

if [ ! -z "${vcf.endsWith(".bcf")?"BCF":""}" ] ; then
	bcftools view -O z -o reformat.vcf.gz "${vcf}"
	bcftools index --tbi reformat.vcf.gz 
	bcftools index --stats reformat.vcf.gz | cut -f 1 > chroms.txt
else

	ln -s "${vcf}" reformat.vcf.gz
	ln -s "${vcf}.tbi" reformat.vcf.gz.tbi
	bcftools index --stats "${vcf}" | cut -f 1 > chroms.txt
fi

##################################################################################
cat <<- EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="vcf">${vcf}</entry>
</properties>
EOF
"""
}


process COMPILE_MINIKIT {
executor "local"
afterScript "rm -rf TMP"
input:
	val(meta)
output:
	path("minikit.jar"),emit:jar
	path("version.xml"),emit:version

script:
"""
${moduleLoad("java/1.8")}

mkdir TMP

cat << "__EOF__" > TMP/Minikit.java
import java.io.BufferedReader;
import java.util.regex.*;
import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.stream.*;

public class Minikit {

private static class MinMax {
		Double min=null;
		Double max= null;
		}

private final Map<String,MinMax> tags = new HashMap<>();

private void fill(String key,String value) {
if(value.equals(".") || value.equalsIgnoreCase("nan")) return;
Double n;
try {
	n = new Double(value);
} catch(Throwable err) {
	return;
}
MinMax mm = tags.get(key);
if(mm.min==null) {
	mm.min= n;
	mm.max= n;
	}
else	{
	if(n.compareTo(mm.min)<0) mm.min=n;
	if(mm.max.compareTo(n)<0) mm.max=n;
	}
}

private void instanceMain(final String args[]) {
	try {
		final Pattern comma = Pattern.compile("[,]");
		final Pattern tab = Pattern.compile("[\t]");
		final Pattern semicolon = Pattern.compile("[;]");
		String mode="";
		int optind=0;
		while(optind < args.length) {
			if(args[optind].equals("--an") && optind+1< args.length) {
				optind++;
				for(String tag : args[optind].split("[, ;]+")) {
					tags.put(tag,new MinMax());
					}
				}
			else if(args[optind].equals("--")) {
				optind++;
				break;
				}
			else if(args[optind].startsWith("-")) {
				System.err.println("unknown option "+args[optind]);
				System.exit(-1);
				}
			optind++;
			}
		if(optind != args.length) {
			System.err.println("illegal number of arguments");
			System.exit(-1);
			}
		try(BufferedReader br = new BufferedReader(new InputStreamReader(System.in))) {
			br.lines().
				filter(L->!L.startsWith("#")).
				map(T->tab.split(T)).
				map(T->T[7]).
				flatMap(T->Arrays.stream(semicolon.split(T))).
				forEach(T->{
					final int i= T.indexOf("=");
					if(i==-1) return;
					final String key = T.substring(0,i);
					if(!tags.containsKey(key)) return;
					final String value = T.substring(i+1);
					if(value.contains(",")) return;
					fill(key,value);
					});
			}
		catch(Exception err) {
			throw err;
			}
		int n=0;
		PrintStream out = System.out;
		for(String key: tags.keySet()) {
			final MinMax mm = tags.get(key);
			if(mm.min==null) continue;
			if(mm.max==null) continue;
			if(mm.min.equals(mm.max)) continue;
			out.print(" -an "+key+" ");
			n++;
			}
		out.println();
		out.flush();
		if(n==0) System.err.println("NO VALID TAGS WAS FOUND IN "+String.join(" , ",tags.keySet()));
	    }
	catch(final Throwable err ) {
		err.printStackTrace();
		System.exit(-1);
		}
	}
public static void main(final String[] args) {
	new Minikit().instanceMain(args);
	}
}
__EOF__


cat << EOF > TMP/tmp.mf
Manifest-Version: 1.0
Main-Class: Minikit
EOF


javac -d TMP -sourcepath TMP TMP/Minikit.java
jar cfm minikit.jar TMP/tmp.mf -C TMP .

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">compile minikit</entry>
	<entry key="javac.version">\$(javac -version 2>&1)</entry>
</properties>
EOF
"""
stub:
"""
touch "minikit.jar"
echo "<properties/>" > version.xml
"""
}


process RECALIBRATE_SNP {
cpus 5
memory "15g"
afterScript 'rm -rf  TMP'
input:
	val(meta)
	val(genomeId)
	val(vcf)
	path(minikit)
output:
	tuple val(vcf),path("snp.recal.vcf.gz"),path("snp.tranches.txt"),emit:output
	path("version.xml"),emit:version
	path("snp.recal.vcf.gz.tbi"),emit:index
	path("snps.plot.R"),emit:R
script:
	def reference = params.genomes[genomeId].fasta
"""
hostname 1>&2
${moduleLoad("gatk/0.0.0 r bcftools")}

mkdir TMP

# 20200630 add AS https://gatk.broadinstitute.org/hc/en-us/articles/360036510892-VariantRecalibrator
# finalement marche pas avec -AS

ANN=\$(bcftools view --type snps -G "${vcf}" | java -jar ${minikit} --an ExcessHet,FS,InbreedingCoeff,QD,ReadPosRankSum,SOR )

test ! -z "\${ANN}"


gatk  --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" VariantRecalibrator  \
	-R "${reference}" \
	-V "${vcf}" \
	\${ANN} \
	${params.vqsr.tranches} \
	--max-gaussians ${params.vqsr.max_gaussians} \
	${params.vqsr[genomeId].recalSnp.trim()} \
	-mode SNP \
	-O "snp.recal.vcf.gz" \
	--tranches-file "snp.tranches.txt" \
        --dont-run-rscript \
	--rscript-file  "snps.plot.R"

##################################################################################
cat <<- EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">VQSR SNP</entry>
	<entry key="vcf">${vcf}</entry>
	<entry key="annotation">\${ANN}</entry>
        <entry key="gatk.version">\$( gatk --version 2> /dev/null  | paste -s -d ' ' )</entry>
</properties>
EOF

"""
}



process RECALIBRATE_INDEL {
memory "15g"
afterScript 'rm -rf  TMP'
cpus 5
input:
	val(meta)
	val(genomeId)
	val(vcf)
	path(minikit)
output:
	tuple val(vcf),path("indel.recal.vcf.gz"),path("indel.tranches.txt"),emit:output
	path("version.xml"),emit:version
	path("indel.recal.vcf.gz.tbi"),emit:index
	path("indels.plot.R"),emit:R
script:
	def reference = params.genomes[genomeId].fasta

"""
hostname 1>&2
${moduleLoad("gatk/0.0.0 r bcftools")}

mkdir -p TMP

ANN=\$(bcftools view --type indels -G "${vcf}" | java -jar ${minikit} --an ReadPosRankSum,InbreedingCoeff,QD,SOR,DP )

test ! -z "\${ANN}"

# 20200630 add AS https://gatk.broadinstitute.org/hc/en-us/articles/360036510892-VariantRecalibrator
# finalement AS ne semble pas marcher avec indel

gatk  --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" VariantRecalibrator  \
	-R "${reference}" \
	-V "${vcf}" \
	\${ANN} \
	${params.vqsr.tranches} \
	--max-gaussians ${params.vqsr.max_gaussians} \
	${params.vqsr[genomeId].recalIndel.trim()} \
	-mode INDEL \
	-O "indel.recal.vcf.gz" \
	--tranches-file "indel.tranches.txt" \
        --dont-run-rscript \
	--rscript-file indels.plot.R

##################################################################################
cat <<- EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">VariantRecalibrator INDEL</entry>
	<entry key="vcf">${vcf}</entry>
	<entry key="annotation">\${ANN}</entry>
        <entry key="gatk.version">\$( gatk --version 2> /dev/null  | paste -s -d ' ' )</entry>
</properties>
EOF

"""
}


process APPLY_RECALIBRATION_SNP {
tag "${row.contig}"
afterScript 'rm -rf  TMP'
memory '10g'
input:
	val(meta)
	val(genomeId)
	val(row)
output:
	tuple val(row),path("recal.vcf.gz"),emit:output
	path("version.xml"),emit:version
	path("recal.vcf.gz.tbi"),emit:index
script:
	def reference = params.genomes[genomeId].fasta

"""
hostname 1>&2
${moduleLoad("gatk/0.0.0 r bcftools")}
mkdir TMP

gatk  --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" ApplyVQSR  \
        -R "${reference}" \
        -V "${row.vcf}"  \
        -mode SNP \
        -L "${row.contig}" \
        --truth-sensitivity-filter-level 99.0 \
        --recal-file "${row.recal_snp_vcf}" \
        --tranches-file "${row.recal_snp_tranches}" \
        -O "TMP/jeter.vcf.gz"

mv TMP/jeter.vcf.gz recal.vcf.gz
mv TMP/jeter.vcf.gz.tbi recal.vcf.gz.tbi

##################################################################################
cat <<- EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">ApplyVQSR SNP</entry>
	<entry key="vcf">${row.vcf}</entry>
	<entry key="contig">${row.contig}</entry>
        <entry key="gatk.version">\$( gatk --version 2> /dev/null  | paste -s -d ' ' )</entry>
</properties>
EOF
"""
}


process APPLY_RECALIBRATION_INDEL {
tag "${row.contig}"
label "gatk"
afterScript 'rm -rf  TMP'
memory '10g'
input:
	val(meta)
	val(genomeId)
	val(row)
output:
	path("recal.vcf.gz"),emit:vcf
	path("recal.vcf.gz.tbi"),emit:index
	path("version.xml"),emit:version
	
script:
	def reference = params.genomes[genomeId].fasta
"""
hostname 1>&2
${moduleLoad("gatk/0.0.0 r bcftools")}
mkdir -p TMP

gatk  --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" ApplyVQSR  \
	-R "${reference}" \
	-V "${row.vcf}"  \
	-mode INDEL \
	-L "${row.contig}" \
	--truth-sensitivity-filter-level 99.0 \
	--recal-file "${row.recal_indel_vcf}" \
	--tranches-file "${row.recal_indel_tranches}" \
	-O "TMP/jeter.vcf.gz"

mv TMP/jeter.vcf.gz recal.vcf.gz
mv TMP/jeter.vcf.gz.tbi recal.vcf.gz.tbi

##################################################################################
cat <<- EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">ApplyVQSR INDEL</entry>
	<entry key="vcf">${row.vcf}</entry>
	<entry key="contig">${row.contig}</entry>
        <entry key="gatk.version">\$( gatk --version 2> /dev/null  | paste -s -d ' ' )</entry>
</properties>
EOF
"""
}
