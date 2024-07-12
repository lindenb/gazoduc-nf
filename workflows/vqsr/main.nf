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
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run my_pipeline --input input_file.csv")
   exit 0
}
// validate parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)


include {moduleLoad;escapeXml} from '../../modules/utils/functions.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include { PICARD_GATHER_VCFS_01 } from '../../modules/picard/picard.gathervcfs.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'




include {BCFTOOLS_CONCAT} from '../../subworkflows/bcftools/bcftools.concat.02.nf'
include {moduleLoad; escapeXml; dumpParams; runOnComplete} from '../../modules/utils/functions.nf'
include {JVARKIT_VCF_TO_INTERVALS_01} from '../../subworkflows/jvarkit/jvarkit.vcf2intervals.nf'

if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}


workflow {
		vcf2bed_ch = BCF_TO_VCF(params.vcf)

		compile_ch= COMPILE_MINIKIT()


		vcf2rgn_ch = JVARKIT_VCF_TO_INTERVALS_01([distance:"10Mb",min_distance:"100"],vcf,output_bed)

		rgn_vcf_ch = vcf2rgn_ch.bed.splitCsv(sep:'\t',header:false).
			map{ [it[0]+":"+((it[1] as int)+1)+"-"+it[2], it[3] ]}


		recal_snp_ch = RECALIBRATE_SNP( genomeId, vcf2bed_ch.output,compile_ch.output, vcf2rgn_ch.vcf2bed, recal_snps)

		recal_indel_ch = RECALIBRATE_INDEL(genomeId, vcf2bed_ch.output ,compile_ch.output, vcf2rgn_ch.vcf2bed, recal_indels)


		

		apply_snp_ch = APPLY_RECALIBRATION_SNP(genomeId , recal_snp_ch.vcf, recal_snp_ch.tranches, rgn_vcf_ch)
		apply_indel_ch = APPLY_RECALIBRATION_INDEL(genomeId , recal_indel_ch.vcf, recal_indel_ch.tranches, apply_snp_ch.output)


		concat_ch = BCFTOOLS_CONCAT([method:"all"],apply_indel_ch.output.map{[vcf:it.toString()]},file("NO_FILE"))		
	}


process BCF_TO_VCF {
tag "${file(vcf).name}"
afterScript 'rm -rf  TMP'
input:
	val(vcf)
output:
	path("reformat.vcf.gz"),emit:output
	path("reformat.vcf.gz.tbi"),emit:index
	path("chroms.txt"),emit:contig
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}
set -o pipefail
mkdir -p TMP

if ${vcf.endsWith(".list")} ; then

	bcftools concat --threads ${task.cpus} --drop-genotypes -a -O z -o TMP/reformat.vcf.gz --file-list "${vcf}"
	bcftools index  --threads ${task.cpus}  --tbi TMP/reformat.vcf.gz 
	mv TMP/reformat.* ./
	
elif ${vcf.endsWith(".bcf")} ; then

	bcftools view  --threads ${task.cpus} --drop-genotypes -O z -o TMP/reformat.vcf.gz "${vcf}"
	bcftools index  --threads ${task.cpus} --tbi TMP/reformat.vcf.gz 
	mv TMP/reformat.* ./

else

	ln -s "${vcf}" reformat.vcf.gz
	ln -s "${vcf}.tbi" reformat.vcf.gz.tbi
fi

bcftools index --stats reformat.vcf.gz | cut -f 1 > chroms.txt

"""
}


process COMPILE_MINIKIT {
executor "local"
afterScript "rm -rf TMP"
output:
	path("minikit.jar"),emit:output
script:
"""
${moduleLoad("java/1.8")}

mkdir -p TMP/TMP
cp -v "${moduleDir}/Minikit.java" TMP/Minikit.java

cat << EOF > TMP/tmp.mf
Manifest-Version: 1.0
Main-Class: Minikit
EOF

javac -d TMP -sourcepath TMP TMP/Minikit.java
jar cfm minikit.jar TMP/tmp.mf -C TMP .
"""
}


process RECALIBRATE_SNP {
tag "file(vcf).name"
afterScript 'rm -rf  TMP'
input:
	val(genomeId)
	val(vcf)
	path(minikit)
	path(bed)
	path(recal_snps)
output:
	path("${params.prefix}snp.recal.vcf.gz"),emit:vcf
	path("${params.prefix}snp.tranches.txt"),emit:tranches
	path("${params.prefix}snp.recal.vcf.gz.tbi")
	path("${params.prefix}snp.tranches.pdf")
script:
	def reference = params.genomes[genomeId].fasta
"""
hostname 1>&2
${moduleLoad("java/1.8 r bcftools")}

mkdir -p TMP

test -s "${recal_snps}"

bcftools view --regions-file "${bed}" --type snps -G "${vcf}" | java -jar ${minikit} --an ExcessHet,FS,InbreedingCoeff,QD,ReadPosRankSum,SOR,MQ,MQRankSum > TMP/args.list
test -s TMP/args.list


cat ${recal_snps} >> TMP/args.list

java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar ${params.gatkjar} \\
	-T VariantRecalibrator \\
	-R "${reference}" \\
	-L "${bed}" \\
	--input "${vcf}" \\
	--arg_file TMP/args.list \\
	-mode SNP \\
	--recal_file "TMP/snp.recal.vcf.gz" \\
	--tranches_file "TMP/snp.tranches.txt" \\
	-rscriptFile  "TMP/snp.plot.R"

mv -v "TMP/snp.recal.vcf.gz" "${params.prefix}snp.recal.vcf.gz"
mv -v "TMP/snp.recal.vcf.gz.tbi" "${params.prefix}snp.recal.vcf.gz.tbi"
mv -v "TMP/snp.tranches.txt" "${params.prefix}snp.tranches.txt"
mv -v "TMP/snp.tranches.txt.pdf" "${params.prefix}snp.tranches.pdf"

"""
}



process RECALIBRATE_INDEL {
tag "${file(vcf).name}"
afterScript 'rm -rf  TMP'
input:
	val(genomeId)
	val(vcf)
	path(minikit)
	path(bed)
	path(recal_indels)
output:
	path("${params.prefix}indel.recal.vcf.gz"),emit:vcf
	path("${params.prefix}indel.tranches.txt"),emit:tranches
	path("${params.prefix}indel.recal.vcf.gz.tbi")
script:
	def reference = params.genomes[genomeId].fasta

"""
hostname 1>&2
${moduleLoad("java/1.8 r bcftools")}
mkdir -p TMP

test -s "${recal_indels}"

bcftools view --regions-file "${bed}" --type indels -G "${vcf}" | java -jar ${minikit} --an ExcessHet,FS,InbreedingCoeff,QD,ReadPosRankSum,SOR,MQ,MQRankSum > TMP/args.list
test -s TMP/args.list


cat ${recal_indels} >> TMP/args.list

java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar ${params.gatkjar} \\
	-T VariantRecalibrator \\
	-R "${reference}" \\
	--input "${vcf}" \\
	-L "${bed}" \\
	--arg_file TMP/args.list \\
	-mode INDEL \\
	--recal_file "TMP/indel.recal.vcf.gz" \
	--tranches_file "TMP/indel.tranches.txt" \
	-rscriptFile  "TMP/indel.plot.R"


mv -v "TMP/indel.recal.vcf.gz" "${params.prefix}indel.recal.vcf.gz"
mv -v "TMP/indel.recal.vcf.gz.tbi" "${params.prefix}indel.recal.vcf.gz.tbi"
mv -v "TMP/indel.tranches.txt" "${params.prefix}indel.tranches.txt"
"""
}


process APPLY_RECALIBRATION_SNP {
tag "${interval} ${vcf}"
afterScript 'rm -rf  TMP'
memory '10g'
input:
	val(genomeId)
	val(recal)
	val(tranches)
	tuple val(interval),val(vcf)
output:
	tuple val(interval),path("recal.snp.vcf.gz"),emit:output
	path("recal.snp.vcf.gz.tbi")
script:
	def reference = params.genomes[genomeId].fasta

"""
hostname 1>&2
${moduleLoad("java/1.8 bcftools")}
mkdir -p TMP

if ${vcf.endsWith(".bcf")} ; then
	bcftools view --regions "${interval}" --threads ${task.cpus} -o TMP/jeter.vcf.gz -O z "${vcf}"
	bcftools index --threads ${task.cpus}  -t TMP/jeter.vcf.gz
fi

java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar ${params.gatkjar} \\
        -T ApplyRecalibration \\
        -R "${reference}" \\
	-L "${interval}" \\
        --input "${vcf.endsWith(".bcf")?"TMP/jeter.vcf.gz":"${vcf}"}" \\
        -mode SNP \\
	--ts_filter_level 99.0 \\
	-recalFile "${recal}" \\
	-tranchesFile "${tranches}" \\
	-o  TMP/recal.snp.vcf.gz

mv -v TMP/recal.* ./

"""
}


process APPLY_RECALIBRATION_INDEL {
tag "${interval} ${vcf}"
afterScript 'rm -rf  TMP'
memory '10g'
input:
	val(genomeId)
	val(recal)
	val(tranches)
	tuple val(interval),val(vcf)
output:
	path("recal.indel.vcf.gz"),emit:output
	path("recal.indel.vcf.gz.tbi")
	
script:
	def reference = params.genomes[genomeId].fasta
"""
hostname 1>&2
${moduleLoad("java/1.8")}
mkdir -p TMP


java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar ${params.gatkjar} \\
        -T ApplyRecalibration \\
        -R "${reference}" \\
	-L "${interval}" \\
        --input "${vcf}" \\
        -mode INDEL \\
	--ts_filter_level 99.0 \\
	-recalFile "${recal}" \\
	-tranchesFile "${tranches}" \\
	-o  TMP/recal.indel.vcf.gz

mv -v TMP/recal.* ./
"""
}

process PLOT_VQSLOD {
afterScript "rm -rf TMP"
tag "${vcf.name}"
input:
	val(meta)
	path(vcf)
output:
	path("${params.prefix}vqsr.pdf"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("r bcftools")} 
mkdir -p TMP

bcftools query -f '%VQSLOD\t%FILTER\t%TYPE\\n' "${vcf}" | sed 's/,OVERLAP\$//' > TMP/jeter.tsv


cat << "__EOF__" > TMP/jeter.R
T1 <- read.table("TMPjeter.tsv",header=FALSE,sep="\t",col.names=c("VQSLOD","FILTER","TYPE"),colClasses=c("numeric","character","character"),stringsAsFactors=FALSE)
head(T1)
pairs = unique(T1[,c("FILTER","TYPE")])
head(pairs)


pdf("${params.prefix}vqsr.pdf")
xminmax = c(-5,5)


for(i in 1:nrow(pairs) ) {
        T2 = T1[T1\$TYPE==pairs[i,"TYPE"] & T1\$FILTER==pairs[i,"FILTER"],]
        d <- density(T2\$VQSLOD)
        plot(d,  xlim = xminmax,
                 ylim = c(0,max(d\$y)),
                xlab = paste(pairs[i,]),
                ylab = "Density",
                main =  paste(pairs[i,]),
                col = ifelse(pairs[i,"FILTER"]=="PASS","green","red")
                ) # plots the results
        }
dev.off()
__EOF__

R --vanilla < TMP/jeter.R

##################################################################################
cat <<- EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">plot VQSL distribution</entry>
	<entry key="vcf">${vcf}</entry>
</properties>
EOF
"""
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

process PLOT_VQSLOD {
afterScript "rm -rf TMP"
tag "${vcf.name}"
input:
	val(meta)
	path(vcf)
output:
	path("${params.prefix}vqsr.pdf"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("r bcftools")} 
mkdir -p TMP

bcftools query -f '%VQSLOD\t%FILTER\t%TYPE\\n' "${vcf}" | sed 's/,OVERLAP\$//' > TMP/jeter.tsv


cat << "__EOF__" > TMP/jeter.R
T1 <- read.table("TMPjeter.tsv",header=FALSE,sep="\t",col.names=c("VQSLOD","FILTER","TYPE"),colClasses=c("numeric","character","character"),stringsAsFactors=FALSE)
head(T1)
pairs = unique(T1[,c("FILTER","TYPE")])
head(pairs)


pdf("${params.prefix}vqsr.pdf")
xminmax = c(-5,5)


for(i in 1:nrow(pairs) ) {
        T2 = T1[T1\$TYPE==pairs[i,"TYPE"] & T1\$FILTER==pairs[i,"FILTER"],]
        d <- density(T2\$VQSLOD)
        plot(d,  xlim = xminmax,
                 ylim = c(0,max(d\$y)),
                xlab = paste(pairs[i,]),
                ylab = "Density",
                main =  paste(pairs[i,]),
                col = ifelse(pairs[i,"FILTER"]=="PASS","green","red")
                ) # plots the results
        }
dev.off()
__EOF__

R --vanilla < TMP/jeter.R

##################################################################################
cat <<- EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">plot VQSL distribution</entry>
	<entry key="vcf">${vcf}</entry>
</properties>
EOF
"""
}
