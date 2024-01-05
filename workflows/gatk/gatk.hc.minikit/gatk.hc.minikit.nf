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
nextflow.enable.dsl=2

def gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults().putReference()

gazoduc.make("bams","NO_FILE").
        description("File containing the paths to the BAM/CRAMS files. One path per line").
	required().
	existingFile().
        put()

gazoduc.make("beds","NO_FILE").
        description("a list of BED files. Variants will be called for each bed file. If not set, the whole genome will be used and split into parts").
        put()

gazoduc.make("dbsnp","").
        description("optional path to dbnsp").
        put()

gazoduc.make("mapq",-1).
        description("mapping quality").
	setInteger().
        put()


gazoduc.make("makewindows_size",100_000).
	description("if no --beds is provided, the reference is split using windows of 'x' bp").
	setInteger().
	put()


gazoduc.make("makewindows_overlap",100).
	description("if no --beds is provided, each window will overlap 'x' bp with the next one").
	setInteger().
	put()


gazoduc.make("gatkjar","/LAB-DATA/BiRD/users/lindenbaum-p/packages/gatk/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar").
        description("path to gatk jar").
	existingFile().
        put()

include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {isBlank;runOnComplete;moduleLoad;getVersionCmd} from '../../../modules/utils/functions.nf'
include {SIMPLE_PUBLISH_01} from '../../../modules/utils/publish.simple.01.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.nf'
include {SCATTER_TO_BED} from '../../../subworkflows/picard/picard.scatter2bed.nf'
include {COLLECT_TO_FILE_01} from '../../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_CONCAT_PER_CONTIG_01} from '../../../subworkflows/bcftools/bcftools.concat.contigs.01.nf'


if( params.help ) {
    gazoduc.usage().
        name("gatk.hc.minikit").
        desc("gatk.hc.minikit").
        print();
    exit 0
} else {
   gazoduc.validate();
}



workflow {
	ch1 = GATK_HC_MINIKIT(params, params.reference, file(params.bams), file(params.beds) )
	html = VERSION_TO_HTML(params,ch1.version)
	}

runOnComplete(workflow);


workflow GATK_HC_MINIKIT {
take:
	meta
	reference
	bams
	beds
main:
	version_ch = Channel.empty()

	if(beds.name.equals("NO_FILE")) {
		scatter_ch = SCATTER_TO_BED(["OUTPUT_TYPE":"ACGT","MAX_TO_MERGE":"1000"], reference) 
		version_ch = version_ch.mix(scatter_ch.version)

		mkwin_ch = MAKE_WINDOWS(meta, scatter_ch.bed)
		version_ch = version_ch.mix(mkwin_ch.version)
		each_bed  = mkwin_ch.output.splitText().map{it.trim()}
		}
	else	{
		each_bed  = Channel.fromPath(beds).splitText().map{it.trim()}
		}


	compile_ch = COMPILE(meta)

	call_ch = PER_BED(meta, reference, compile_ch.jar, bams, each_bed)
	version_ch = version_ch.mix(call_ch.version)

	ch1 = COLLECT_TO_FILE_01([:], call_ch.output.collect())
	version_ch = version_ch.mix(ch1.version)

	ch2 = BCFTOOLS_CONCAT_PER_CONTIG_01([:], ch1.output)
	version_ch = version_ch.mix(ch2.version)

        version_ch = MERGE_VERSION(meta, "hc-gatk-minikit", "hc-gatk-minikit",version_ch.collect())

emit:
        vcfs = ch2.vcfs
        version= version_ch

}


process MAKE_WINDOWS {
tag "${bed}"
executor "local"
input:
	val(meta)
	path(bed)
output:
	path("windows.beds.list"),emit:output
	path("version.xml"),emit:version
script:
	def w = ((meta.makewindows_size?:100_000) as int)
	def s = w - ((meta.makewindows_overlap?:100) as int)
"""
hostname 1>&2
${moduleLoad("bedtools")}
set -o pipefail

mkdir -p BEDS

bedtools makewindows -b "${bed}" -w ${w} -s ${s} |\
	grep -vE '^(hs37d5|chrM|chrMT|MT|chrEBV)' |\
	LC_ALL=C sort -T . -t '\t' -k1,1V -k2,2n |\
	split -a 9  --lines=1 --additional-suffix=.bed - BEDS/window

find \${PWD}/BEDS -type f -name "window*.bed"  > windows.beds.list

test -s windows.bed.list

sleep 5

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">split bed genome into parts</entry>
        <entry key="bed">${bed}</entry>
        <entry key="win.size">${w}</entry>
        <entry key="shift.size">${s}</entry>
	<entry key="versions">${getVersionCmd("bedtools")}</entry>
</properties>
EOF
"""
}

process COMPILE {
executor "local"
input:
	val(meta)
output:
	path("minikit.jar"),emit:jar
	path("version.xml"),emit:version
script:
"""
${moduleLoad("openjdk/11.0.8")}
mkdir -p TMP

cat << __EOF__ > Minikit.java


import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.samtools.util.IOUtil;

import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.*;


public class Minikit extends Main{
    private static final Logger LOG = LogManager.getLogger(Main.class);
    @Argument(fullName = StandardArgumentDefinitions.INTERVALS_LONG_NAME, shortName = StandardArgumentDefinitions.INTERVALS_SHORT_NAME, doc = "interval as BED file", common = true, optional = false)
	private String bedPath = null;
    @Argument(fullName = StandardArgumentDefinitions.REFERENCE_LONG_NAME, shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME, doc = "Indexed fasta reference", common = true, optional = false)
	private String reference=null;
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "File to which variants should be written")
	private String finalVcf = null;
    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, doc = "A file containing the path to the bams", common = true, optional = false)
    private String bamsList;
    @Argument(fullName = "minimum-mapq", shortName = "Q", doc = "minimum mapping quality", common = true, optional = true)
	private int mapq= 10;
    @Argument(fullName = "dbsnp", shortName = "dbsnp", doc = "dbsnp", common = true, optional = true)
	private String dbsnp = null;

	private void execute(final List<String> argv) throws Exception {
		final String[] args = argv.toArray(new String[argv.size()]);
		LOG.info(getCommandLineName()+":Executing:  gatk "+ String.join(" ",argv));
		final CommandLineProgram program =
			this.setupConfigAndExtractProgram(args, 
				this.getPackageList(),
				this.getClassList(),
				this.getCommandLineName()
				);
	    final Object result = Main.runCommandLineProgram(program, args);
		if(result==null) return;
		if(Boolean.TRUE.equals(result)) return;
		LOG.warn("Returned "+ result.getClass());
		LOG.error("Result is "+ result);
		final Throwable err= (result instanceof Throwable?Throwable.class.cast(result):null);
		if(err!=null) {
			throw new RuntimeException(err);
			}
		else
			{
			throw new RuntimeException("Failure");
			}
		}


	private void fill_cmd(List<String> cmd,final Path tmpDir,final String bed) {
		cmd.add("-R");
		cmd.add(this.reference.toString());
		//cmd.add("--verbosity");
		//cmd.add("ERROR");
		cmd.add("-L");
		cmd.add(bed);
		cmd.add("--tmp-dir");
		cmd.add(tmpDir.toString());
		if(!StringUtil.isBlank(this.dbsnp)) {
			cmd.add("--dbsnp");
			cmd.add(String.valueOf(this.dbsnp));
			}
		}

	private Path invoke_gvcf(final Path tmpDir,int idx,final String bam,final String bed) throws Exception {
		final Path outPath = tmpDir.resolve(String.format("hc.%03d.g" + FileExtensions.COMPRESSED_VCF, idx));
		final List<String> cmd= new ArrayList<>();
		cmd.add("HaplotypeCaller");
		fill_cmd(cmd,tmpDir,bed);
		cmd.add("-I");
		cmd.add(bam);
		cmd.add("-O");
		cmd.add(outPath.toString());
		cmd.add("-ERC");
		cmd.add("GVCF");
		if(this.mapq>0) {
			cmd.add("--minimum-mapping-quality");
			cmd.add(String.valueOf(this.mapq));
			}
		
		cmd.add("-G");cmd.add("StandardAnnotation");
		cmd.add("-G");cmd.add("AS_StandardAnnotation");
		cmd.add("-G");cmd.add("StandardHCAnnotation");
		
		execute(cmd);
		return outPath;
		}

	private Path invoke_combine(final Path tmpDir,int idx,final List<Path> gvcfs,final String bed) throws Exception {
		final Path outPath = tmpDir.resolve(String.format("combine.%03d.g" + FileExtensions.COMPRESSED_VCF, idx));
		final List<String> cmd= new ArrayList<>();
		cmd.add("CombineGVCFs");
		fill_cmd(cmd,tmpDir,bed);
		for(Path gvcf:gvcfs) {
			cmd.add("-V");
			cmd.add(gvcf.toString());
			}
		cmd.add("-O");
		cmd.add(outPath.toString());
		
		cmd.add("-G");cmd.add("StandardAnnotation");
		cmd.add("-G");cmd.add("AS_StandardAnnotation");

		
		execute(cmd);
		return outPath;
		}
	
	private void invoke_genotype(final Path tmpDir,final Path gvcf,final String bed,final String outPath) throws Exception {
		final List<String> cmd= new ArrayList<>();
		cmd.add("GenotypeGVCFs");
		fill_cmd(cmd,tmpDir,bed);
		cmd.add("-V");
		cmd.add(gvcf.toString());
		cmd.add("-G");cmd.add("StandardAnnotation");
		cmd.add("-G");cmd.add("AS_StandardAnnotation");
		cmd.add("-O");
		cmd.add(outPath);
		execute(cmd);
		}


private int doWork(final String[] args) {
	try {
		final  CommandLineArgumentParser cmdLineParser  = new CommandLineArgumentParser(this);
        final boolean ret = cmdLineParser.parseArguments(System.err, args);
        if (!ret) {
        	System.err.println(cmdLineParser.usage(false, false));
        	return -1;
        	}

		IOUtil.assertFileIsReadable(Paths.get(this.bedPath));
		
		//final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(Paths.get(this.reference));
		//final GenomeLocParser parser = new GenomeLocParser(dict);
		//final Locatable region  = parser.parseGenomeLoc(this.bedPath);
		
		List<String> bams = Files.readAllLines(Paths.get(this.bamsList));
		if (bams.isEmpty()) {
        	LOG.error("No Bam was provided");
        	return -1;
        	}
		
		
		Path tmpDir = Files.createTempDirectory("tmp");
		final List<Path> gvcfs_list = new ArrayList<>(bams.size());
		for(int i=0;i< bams.size();i++) {
			LOG.info("("+(i+1)+"/"+bams.size()+" "+bams.get(i));
			final Path g_vcf_gz = invoke_gvcf(tmpDir,i,bams.get(i),this.bedPath);
			gvcfs_list.add(g_vcf_gz);
			}
			
		final int sqrt = Math.max(20, (int)Math.sqrt(gvcfs_list.size()));
		
		final List<Path> combined0_list = new ArrayList<>(sqrt);
		int i=0;
		while(!gvcfs_list.isEmpty()) {
			final List<Path> L = new ArrayList<>(sqrt);
			while(!gvcfs_list.isEmpty() && L.size() < sqrt) {
				L.add(gvcfs_list.remove(0));
				}
			if(L.size()==1) {
				combined0_list.add(L.get(0));
				}
			else {
				final Path combined = invoke_combine(tmpDir,i,L, this.bedPath);
				combined0_list.add(combined);
				}
			i++;
			}
		Path to_genotype;
		if(combined0_list.size()==1) {
			to_genotype = combined0_list.get(0);
			}
		else
			{
			to_genotype = invoke_combine(tmpDir,i,combined0_list,this.bedPath);
			i++;
			}
		invoke_genotype(tmpDir,to_genotype,bedPath,this.finalVcf);
		return 0;
		}
	catch(Throwable err) {
		err.printStackTrace();
		return -1;
		}
	}

@Override
protected String getCommandLineName() {
	return this.getClass().getSimpleName();
	}


public static void main(String[] args)
    {
	int ret= new Minikit().doWork(args);
	System.exit(ret);
    }
}


__EOF__


javac -d TMP -cp ${meta.gatkjar} -sourcepath . Minikit.java
jar cvf minikit.jar -C TMP .
rm -rf TMP

sleep 5

#####
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">compile gatk</entry>
</properties>
EOF
"""
}


process PER_BED {
afterScript "rm -rf TMP"
tag "${file(bed).name}"
memory "10g"
cpus 3
input:
	val(meta)
	val(reference)
	path(minikit)
	path(bams)
	val(bed)
output:
	path("genotyped.bcf"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("openjdk/11.0.8 bcftools/0.0.0")}

mkdir -p TMP

java -Xmx${task.memory.giga}g -Dsamjdk.compression_level=1 -Djava.io.tmpdir=TMP -cp ${meta.gatkjar}:${minikit} Minikit \
                -I "${bams}" \
                -L "${bed}" \
                ${meta.mapq && ((meta.mapq as int)>0)?"--minimum-mapq ${meta.mapq}":""} \
                --output TMP/selection.vcf.gz \
                --reference "${reference}" \
                ${isBlank(meta.dbsnp)?"":"--dbsnp ${meta.dbsnp}"}


bcftools view --threads ${task.cpus} --compression-level 9 -O b -o genotyped.bcf TMP/selection.vcf.gz 
bcftools index -f --threads ${task.cpus} genotyped.bcf

#####
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">call bed</entry>
        <entry key="bed">${bed}</entry>
        <entry key="bams">${bams}</entry>
        <entry key="mapq">${meta.mapq?:""}</entry>
        <entry key="dbsnp">${meta.dbsnp}</entry>
</properties>
EOF
"""
}
