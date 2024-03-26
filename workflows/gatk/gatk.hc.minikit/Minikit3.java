
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;


public class Minikit {
	private static final String COMPRESSED_VCF = ".vcf.gz";
	
    private final static Logger LOG = Logger.getLogger(Minikit.class);
	private String bedPath = null;
	private String reference=null;
	private String finalVcf = null;
    private String bamsList;
	private int mapq= 10;
	private String dbsnp = null;
	private int attempt=1;
	private int cpus=1;

	private void execute(final List<String> argv) throws Exception {
	    final org.broadinstitute.gatk.engine.CommandLineGATK instance= new org.broadinstitute.gatk.engine.CommandLineGATK();
		final String[] args = argv.toArray(new String[argv.size()]);
		LOG.info(":Executing:  gatk "+ String.join(" ",argv));
		
		org.broadinstitute.gatk.engine.CommandLineGATK.start(instance, args);
		
	    final int result =org.broadinstitute.gatk.engine.CommandLineGATK.result;
		LOG.info("Result is "+ result);
		System.gc();
		if(result!=0) {
			throw new RuntimeException("Failure "+result);
			}
		}


	private void fill_cmd(List<String> cmd,final Path tmpDir,final String bed) {
		cmd.add("-R");
		cmd.add(this.reference.toString());
		//cmd.add("--verbosity");
		//cmd.add("ERROR");
		cmd.add("-L");
		cmd.add(bed);
		if(!StringUtil.isBlank(this.dbsnp)) {
			cmd.add("--dbsnp");
			cmd.add(String.valueOf(this.dbsnp));
			}
		}

	private Path invoke_gvcf(final Path tmpDir,int idx,final String bam,final String bed) throws Exception {
		final Path outPath = tmpDir.resolve(String.format("hc.%03d.g" + COMPRESSED_VCF, idx));
		final List<String> cmd= new ArrayList<>();
		cmd.add("-T");
		cmd.add("HaplotypeCaller");
		fill_cmd(cmd,tmpDir,bed);
		cmd.add("-I");
		cmd.add(bam);
		cmd.add("-o");
		cmd.add(outPath.toString());
		cmd.add("-ERC");
		cmd.add("GVCF");
		if(this.mapq>0) {
			cmd.add("--min_mapping_quality_score");
			cmd.add(String.valueOf(this.mapq));
			}
		if(this.cpus>1) {
			cmd.add("-nct");
			cmd.add(String.valueOf(this.cpus));
			}
		
		cmd.add("-G");cmd.add("StandardAnnotation");
		cmd.add("-G");cmd.add("AS_StandardAnnotation");
		cmd.add("-G");cmd.add("StandardHCAnnotation");
		
		execute(cmd);
		return outPath;
		}

	private Path invoke_combine(final Path tmpDir,int idx,final List<Path> gvcfs,final String bed) throws Exception {
		final Path outPath = tmpDir.resolve(String.format("combine.%03d.g" +  COMPRESSED_VCF, idx));
		final List<String> cmd= new ArrayList<>();
		cmd.add("-T");
		cmd.add("CombineGVCFs");
		fill_cmd(cmd,tmpDir,bed);
		for(Path gvcf:gvcfs) {
			cmd.add("-V");
			cmd.add(gvcf.toString());
			}
		cmd.add("-o");
		cmd.add(outPath.toString());
		
		cmd.add("-G");cmd.add("StandardAnnotation");
		cmd.add("-G");cmd.add("AS_StandardAnnotation");

		
		execute(cmd);
		return outPath;
		}
	
	private void invoke_genotype(final Path tmpDir,final Path gvcf,final String bed,final String outPath) throws Exception {
		final List<String> cmd= new ArrayList<>();
		cmd.add("-T");
		cmd.add("GenotypeGVCFs");
		fill_cmd(cmd,tmpDir,bed);
		cmd.add("-V");
		cmd.add(gvcf.toString());
		cmd.add("-G");cmd.add("StandardAnnotation");
		cmd.add("-G");cmd.add("AS_StandardAnnotation");
		cmd.add("-o");
		cmd.add(outPath);
		execute(cmd);
		}


private int doWork(final String[] args) {
	try {
	    LOG.setLevel(org.apache.log4j.Level.ALL);
        int optind=0;
        while(optind<args.length) {
        	String arg = args[optind];
        	optind++;
        	if(arg.equals("-L")) {
        		this.bedPath = args[optind++];
        		}
        	else if(arg.equals("--dbsnp")) {
        		this.dbsnp = args[optind++];
        		}
        	else if(arg.equals("-R")) {
        		this.reference = args[optind++];
        		}
        	else if(arg.equals("-I")) {
        		this.bamsList = args[optind++];
        		}
        	else if(arg.equals("--mapq")) {
        		this.mapq = Integer.parseInt(args[optind++]);
        		}
        	else if(arg.equals("-o")) {
        		this.finalVcf = args[optind++];
        		}
        	else if(arg.equals("--attempt")) {
        		this.attempt = Integer.parseInt(args[optind++]);
        		}
        	else if(arg.equals("--threads")) {
        		this.cpus = Integer.parseInt(args[optind++]);
        		}
        	else
        		{
        		LOG.error("illegal number of args ("+arg+")");
        		return -1;
        		}
        	}
		if(bamsList==null) {
			LOG.error("bam missing");
			return -1;
			}
		if(bedPath==null) {
			LOG.error("bam missing");
			return -1;
			}
		if(finalVcf==null) {
			LOG.error("output missing");
			return -1;
			}
		if(reference==null) {
			LOG.error("output missing");
			return -1;
			}

		IOUtil.assertFileIsReadable(Paths.get(this.bedPath));
		
		//final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(Paths.get(this.reference));
		//final GenomeLocParser parser = new GenomeLocParser(dict);
		//final Locatable region  = parser.parseGenomeLoc(this.bedPath);
		
		final List<String> bams = Files.readAllLines(Paths.get(this.bamsList));
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


public static void main(String[] args)
    {
	int ret= new Minikit().doWork(args);
	System.exit(ret);
    }
}
