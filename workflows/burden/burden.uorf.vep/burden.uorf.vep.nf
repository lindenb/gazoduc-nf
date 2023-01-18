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
nextflow.enable.dsl=2

include {getGnomadGenomePath;isHg19;runOnComplete;moduleLoad;getKeyValue;hasFeature;getVersionCmd} from '../../../modules/utils/functions.nf'
include {UTR_ANNOTATOR_DOWNLOAD_01} from '../../../modules/vep/utrannotator.download.01.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.nf'
include {RVTESTS_REHEADER_01} from '../../../modules/rvtests/rvtests.reheader.01.nf'
include {RVTESTS_POST_PROCESS} from '../../../subworkflows/rvtests/rvtests.post.process.01.nf'
include {CONCAT_FILES_01} from '../../../modules/utils/concat.files.nf'
include {VCF_TO_BED} from '../../../modules/bcftools/vcf2bed.01.nf'
include {VALIDATE_CASE_CONTROL_PED_01} from '../../../modules/pedigree/validate.case.ctrl.pedigree.01.nf'
include {VCF_INTER_CASES_CONTROLS_01} from '../../../subworkflows/bcftools/vcf.inter.cases.controls.01.nf'
include {PEDIGREE_FOR_RVTESTS} from '../../../modules/rvtests/rvtests.cases.controls.ped.01.nf'
include {RVTESTS01_VCF_01} from '../../../modules/rvtests/rvtests.vcf.01.nf'

params.reference=""
params.pedigree=""
params.vcf=""
params.disableFeatures="";
params.help=false
params.kozak_strength="Weak,Moderate,Strong"
params.gnomadPop="AF_nfe"
params.gnomadAF = 0.01


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
nextflow -C ../../confs/cluster.cfg  run -resume burden.first.intron.01.nf \\
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
		burden_ch = BURDEN_UPSTREAM_ORF_VEP(params, params.reference, params.vcf, file(params.pedigree))
		ZIPIT(params,burden_ch.zip.collect())
		}

workflow BURDEN_UPSTREAM_ORF_VEP {
	take:
		meta
		reference
		vcf
		pedigree
	main:

		version_ch = Channel.empty()
		to_zip = Channel.empty()
		
		vcfbed_ch = VCF_TO_BED(meta,vcf)
		version_ch = version_ch.mix(vcfbed_ch.version)

		plugin_ch = UTR_ANNOTATOR_DOWNLOAD_01(meta,reference)
		version_ch = version_ch.mix(plugin_ch.version)
		
		intersect_ch = INTERSECT(meta, vcfbed_ch.bed, plugin_ch.bed)
		version_ch = version_ch.mix(intersect_ch.version)
		
		minikit_ch = COMPILE_MINIKIT(meta)
		version_ch = version_ch.mix(minikit_ch.version)

		ped_ch = VALIDATE_CASE_CONTROL_PED_01(meta,pedigree)
		version_ch = version_ch.mix(ped_ch.version)
		
		rvtests_ped_ch = PEDIGREE_FOR_RVTESTS(meta,ped_ch.pedigree)
		version_ch = version_ch.mix(rvtests_ped_ch.version)

		each_bed = intersect_ch.bed.splitCsv(header:false,sep:'\t').map{T->[T[0],T[3],T[4]]}

		annot_ch = ANNOTATE(meta, reference, minikit_ch.jar, rvtests_ped_ch.pedigree, plugin_ch.output, each_bed)
		version_ch = version_ch.mix(annot_ch.version)
		
		assoc_ch = RVTESTS01_VCF_01(meta, reference, annot_ch.vcf, rvtests_ped_ch.pedigree)
		version_ch = version_ch.mix(assoc_ch.version)

		concat_ch = CONCAT_FILES_01(meta,assoc_ch.output.collect())
		version_ch = version_ch.mix(concat_ch.version)

		digest_ch = RVTESTS_POST_PROCESS(meta, reference,file(vcf).name,concat_ch.output)
                version_ch = version_ch.mix(digest_ch.version)
		to_zip = to_zip.mix(digest_ch.zip)

		version_ch = MERGE_VERSION(meta, "BurdenMicroORF", "Burden in micro-ORF.", version_ch.collect())
		to_zip = to_zip.mix(version_ch.version)

		html = VERSION_TO_HTML(params,version_ch.version)
		to_zip = to_zip.mix(html.html)
		
	emit:
		version = version_ch
		zip = to_zip
	}

process INTERSECT {
executor "local"
input:
	val(meta)
	path(vcfbed)
	path(pluginbed)
output:
	path("intersect.bed"),emit:bed
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bedtools")}

sort -T . -t '\t' -k1,1 -k2,2n "${vcfbed}" > jeter.a.bed
sort -T . -t '\t' -k1,1 -k2,2n "${pluginbed}" > jeter.b.bed

bedtools intersect -u -a jeter.a.bed -b jeter.b.bed |\
	awk -F '\t' '{printf("%s\t${pluginbed.toRealPath()}\\n",\$0);}' > intersect.bed

rm -f jeter.a.bed jeter.b.bed


###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">intersection vcf.bed/utr.bed</entry>
	<entry key="versions">${getVersionCmd("bedtools")}</entry>
</properties>
EOF
"""
}

process COMPILE_MINIKIT {
	executor "local"
	input:
		val(meta)
	output:
		path("minikit.jar"),emit:jar
		path("version.xml"),emit:version
	script:
		def strengths = meta.kozak_strength?:"Weak,Moderate,Strong"
	"""
	hostname 1>&2
	${moduleLoad("picard/0.0.0 java-jdk/8.0.112")}
	mkdir -p TMP

cat << __EOF__ > TMP/MiniKit.java
import java.io.*;
import java.nio.file.*;
import java.nio.charset.*;
import java.util.*;
import java.util.regex.*;
import java.util.stream.*;
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.*;
import htsjdk.variant.vcf.*;


public class MiniKit {

void instanceMainWithExit(final String args[]) {
	try {
		final Pattern pipe = Pattern.compile(Pattern.quote("|"));
		final Pattern comma = Pattern.compile(Pattern.quote(","));
		final Pattern amp = Pattern.compile(Pattern.quote("&"));
		final Set<String> kz_strengths = Arrays.stream(comma.split("${strengths}")).collect(Collectors.toSet());
		try(final VCFIterator in = new VCFIteratorBuilder().open(System.in)) {

		final VCFHeader h = in.getHeader();
		final VCFInfoHeaderLine csqInfo = h.getInfoHeaderLine("CSQ");
		if(csqInfo==null) throw new IllegalArgumentException("CSQ not found");
		String description = csqInfo.getDescription();
		int colon = description.indexOf(":");
		description= description.substring(colon+1);
		final String[] cols = pipe.split(description);
		int feature_column = -1;
		int csq_column = -1;
		for(int i=0;i< cols.length;i++) {
			cols[i]=cols[i].trim();
			if(cols[i].equals("Feature")) feature_column = i ;
			else if(cols[i].equals("five_prime_UTR_variant_annotation")) csq_column = i;
			}
	
		if(feature_column==-1) throw new IllegalArgumentException("Cannot find CSQ/Feature in "+ description);
		if(csq_column==-1) throw new IllegalArgumentException("Cannot find CSQ/five_prime_UTR_variant_annotation in "+description);

		final VariantContextWriterBuilder vcwb=new VariantContextWriterBuilder();
		vcwb.setCreateMD5(false);
		vcwb.setOutputStream(System.out);
		vcwb.setReferenceDictionary(null);
		vcwb.clearOptions();
		final VariantContextWriter w = vcwb.build();
		
		w.writeHeader(h);
		while(in.hasNext()) {
			final VariantContext ctx = in.next();
			final List<String> annot_in = ctx.getAttributeAsStringList("CSQ","");
			final List<String> annot_out = new ArrayList<>(annot_in.size());
			for(final String ann: annot_in) {
				final String[] tokens= pipe.split(ann);
				if(tokens.length<=feature_column) continue;
				if(tokens.length<=csq_column) continue;
				if(tokens[feature_column].isEmpty()) continue;
				if(tokens[csq_column].isEmpty()) continue;
				final Map<String,String> hash = new HashMap<>();
				
				for(final String s : amp.split(tokens[csq_column]) ) {
					colon = s.indexOf(":");
					if(colon==-1) continue;
					final String key = s.substring(0,colon);
					final String value = s.substring(colon+1);
					hash.put(key,value);
					}
				if(hash.entrySet().
					stream().
					filter(KV->KV.getKey().endsWith("KozakStrength")).
					map(KV->KV.getValue()).
					noneMatch(V->kz_strengths.contains(V))
					) continue;

				annot_out.add(ann);
				}
			if(annot_out.isEmpty()) continue;
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			vcb.attribute("CSQ",annot_out);
			w.add(vcb.make());
			}
		w.close();
		}//end iterator
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

javac -cp "\${PICARD_JAR}" -d TMP -sourcepath TMP TMP/MiniKit.java
echo "Manifest-Version: 1.0" > TMP/tmp.mf
echo "Main-Class: MiniKit" >> TMP/tmp.mf
echo "Class-Path: \${PICARD_JAR}" | fold -w 71 | awk '{printf("%s%s\\n",(NR==1?"": " "),\$0);}' >>  TMP/tmp.mf
jar cvfm minikit.jar TMP/tmp.mf -C TMP .


###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">compile minikit</entry>
	<entry key="kozak.strengths">${strengths}</entry>
	<entry key="picard.jar">\${PICARD_JAR}</entry>
	<entry key="versions">${getVersionCmd("javac")}</entry>
</properties>
EOF



"""
}




process ANNOTATE {
tag "${contig} ${vcf}"
afterScript "rm -rf TMP"
memory "5g"
input:
	val(meta)
	val(reference)
	path(minikit)
	path(rvtestpedigree)
	path(pluginfile)
	tuple val(contig),val(vcf),val(bedall)
output:
	path("contig.bcf"),emit:vcf
	path("version.xml"),emit:version
script:
	def vep_module=(isHg19(reference)?"ensembl-vep/104.3":"")
	def gnomadgenome= getGnomadGenomePath(meta,reference)
	def gnomadgenomefilterexpr = "FILTER~\"GNOMAD_GENOME_BAD_AF\"|| FILTER~\"GNOMAD_GENOME_InbreedingCoeff\"|| FILTER~\"GNOMAD_GENOME_RF\""
	def gnomadPop = meta.gnomadPop?:"AF_nfe"
	def gnomadAF = meta.gnomadAF?:0.01
"""
hostname 1>&2
${moduleLoad("bcftools jvarkit "+vep_module)}
set -x
mkdir -p TMP

awk -F '\t' '(\$1=="${contig}")' "${bedall}" > TMP/jeter.bed


# case/list for contrast
awk -F '\t' '(\$6=="1") {printf("%s\\n",\$2);}' "${rvtestpedigree}" > TMP/ctrl.list
awk -F '\t' '(\$6=="2") {printf("%s\\n",\$2);}' "${rvtestpedigree}" > TMP/cases.list
cat TMP/ctrl.list TMP/cases.list > TMP/samples.txt

bcftools view -m 2 -M 3 --samples-file TMP/samples.txt  -O u --regions-file TMP/jeter.bed "${vcf}" -o TMP/jeter2.bcf
mv TMP/jeter2.bcf TMP/jeter1.bcf

bcftools view -m 2 -M 3 --min-ac 1 -O u -o TMP/jeter2.bcf TMP/jeter1.bcf
mv TMP/jeter2.bcf TMP/jeter1.bcf

bcftools norm -f "${reference}" --multiallelics -any -O u -o TMP/jeter2.bcf TMP/jeter1.bcf
mv TMP/jeter2.bcf TMP/jeter1.bcf

set +u
export PERL5LIB=\${PERL5LIB}:${file(pluginfile.toRealPath()).getParent()} 

bcftools view TMP/jeter1.bcf |\
	vep --plugin UTRannotator,${pluginfile.toRealPath()} --cache --format vcf --force_overwrite \
		--output_file STDOUT \
		--no_stats --offline  --dir_cache /LAB-DATA/BiRD/resources/apps/vep  \
		--species homo_sapiens --cache_version 91 \
		--assembly GRCh37  --fasta "${reference}" --use_given_ref --vcf |\
	java -jar "${minikit}" |\
	bcftools view -O u -o TMP/jeter2.bcf -
mv TMP/jeter2.bcf TMP/jeter1.bcf

bcftools +contrast -0 TMP/ctrl.list -1 TMP/cases.list -a PASSOC,FASSOC,NASSOC,NOVELAL,NOVELGT -O u -o TMP/jeter2.bcf TMP/jeter1.bcf
mv TMP/jeter2.bcf TMP/jeter1.bcf

# gnomad genome
bcftools view TMP/jeter1.bcf |\
	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/vcfgnomad.jar --bufferSize 1000 -F '${gnomadPop}' -g "${gnomadgenome}" --max-af "${gnomadAF}" |\
	bcftools view -e '${gnomadgenomefilterexpr}' -O u -o TMP/jeter2.bcf
mv TMP/jeter2.bcf TMP/jeter1.bcf


bcftools sort -T TMP -O b -o contig.bcf TMP/jeter1.bcf
bcftools index contig.bcf

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">annotate variant with vep+UTRannotator</entry>
	<entry key="gnomad.vcf">${gnomadgenome}</entry>
	<entry key="gnomad.pop">${gnomadPop}</entry>
	<entry key="gnomad.AF">${gnomadAF}</entry>
	<entry key="versions">${getVersionCmd("bcftools awk vep jvarkit/vcfgnomad")}</entry>
</properties>
EOF
"""
}



process RVTEST_CONTIG {
tag "${bed.name}"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
	path(vcfs)
	path(pedigree)
	path(reheader)
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

cut -f 1,2,3 "${bed}" | sort -T TMP -t '\t' -k1,1 -k2,2n | bedtools merge > TMP/jeter.bed
bcftools concat -a --regions-file TMP/jeter.bed --file-list "${vcfs}" -O b -o TMP/jeter.bcf
bcftools index TMP/jeter.bcf

i=1
awk -F '\t' '{printf("%s:%d-%s\t%s\\n",\$1,int(\$2)+1,\$3,\$4);}' "${bed}" | while read RGN ENST
do
	bcftools view TMP/jeter.bcf "\${RGN}" |\
		java -jar \${JVARKIT_DIST}/vcfburdenfiltergenes.jar -a "\${ENST}" |\
		bcftools annotate -O z --rename-chrs "${reheader}" -o TMP/jeter.vcf.gz -

	if test `bcftools query -f '.' TMP/jeter.vcf.gz | wc -c` -gt 0 ; then

		bcftools index --tbi -f TMP/jeter.vcf.gz
	

		echo -n "## \${ENST}: " && bcftools query -f '.' TMP/jeter.vcf.gz | wc -c

		# build setFile
		echo "\${ENST}\t\${RGN}" | sed 's/\tchr/\t/' > TMP/variants.setfile

		rvtest  --noweb \
        		--inVcf TMP/jeter.vcf.gz \
			--setFile TMP/variants.setfile \
	        	--pheno "${pedigree}" \
		        --out "ASSOC/part.\${i}" \
			${rvtest_params} 1>&2 2> TMP/last.rvtest.log

		i=\$((i+1))

	fi
done

find \${PWD}/ASSOC -type f -name "part.*assoc" > assoc.list

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">invoke rvtest for bed</entry>
	<entry key="rvtest.path">\$(which rvtest)</entry>
	<entry key="versions">${getVersionCmd("bedtools rvtest bcftools jvarkit/vcfburdenfiltergenes")}</entry>
	<entry key="bed">${bed}</entry>
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
