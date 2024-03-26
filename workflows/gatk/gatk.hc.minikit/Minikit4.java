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


