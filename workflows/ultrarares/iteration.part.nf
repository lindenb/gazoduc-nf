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
include {GATK_BAM2VCF                      } from '../../modules/gatk/bam2vcf'
include {JVARKIT_VCFGNOMAD                 } from '../../modules/jvarkit/vcfgnomad'
include {JVARKIT_VCFFILTERJDK              } from '../../modules/jvarkit/vcffilterjdk'
include {SNPEFF_APPLY                      } from '../../modules/snpeff/apply'
include { BCFTOOLS_SORT                    } from '../../modules/bcftools/sort'
include { BCFTOOLS_INDEX                   } from '../../modules/bcftools/index'
include { BCFTOOLS_CONCAT                  } from '../../modules/bcftools/concat3'
include { VCF_TO_BED                       } from '../../modules/bcftools/vcf2bed2'
include { BED_CLUSTER                      } from '../../modules/jvarkit/bedcluster'
include { BEDTOOLS_SLOP                    } from '../../modules/bedtools/slop'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE1} from '../../modules/bedtools/merge'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE2} from '../../modules/bedtools/merge'



workflow ULTRA_RARES_ITERATION {
take:
	metadata
	fasta
	fai
	dict
	snpeff_db // [meta,directory,name]
	pedigree_ch // [meta,pedigree]
	gnomad_vcf //[meta,vcf,tbi]
	jvarkit_filter // [meta,filter]
	bed_in
	bam_files
main:
	versions = Channel.empty()
	multiqc = Channel.empty()

	GATK_BAM2VCF(
		fasta,
		fai,
		dict,
		[[id:"no_dbsnp"],[],[]],
		pedigree_ch,
		[[id:"no_extra_ref"],[]],
		bam_files.combine(bed_in)
			.map{meta1,hts_files,_meta2,bed->[meta1.plus(id:bed.baseName),hts_files,bed]}
		)
	versions = versions.mix(GATK_BAM2VCF.out.versions)
	vcfs = GATK_BAM2VCF.out.vcf.map{meta,vcf,tbi,bed->[meta,vcf]}
	
	JVARKIT_VCFFILTERJDK(
		jvarkit_filter,
		[[id:"no_ped"],[]],
		vcfs
		)
	versions = versions.mix(JVARKIT_VCFFILTERJDK.out.versions)
	vcfs = JVARKIT_VCFFILTERJDK.out.vcf

	JVARKIT_VCFGNOMAD(
		gnomad_vcf,
		vcfs
		)
	versions = versions.mix(JVARKIT_VCFGNOMAD.out.versions)
	vcfs = 	JVARKIT_VCFGNOMAD.out.vcf
	
	SNPEFF_APPLY(
		snpeff_db,
		vcfs
		)
	versions = versions.mix(SNPEFF_APPLY.out.versions)
	vcfs = 	SNPEFF_APPLY.out.vcf
	
	BCFTOOLS_SORT(vcfs)
	versions = versions.mix(BCFTOOLS_SORT.out.versions)
	vcfs = 	BCFTOOLS_SORT.out.vcf

	BCFTOOLS_INDEX(vcfs)
	versions = versions.mix(BCFTOOLS_INDEX.out.versions)
	vcfs = 	BCFTOOLS_INDEX.out.vcf

	BCFTOOLS_CONCAT(
		vcfs.flatMap{meta,vcf,tbi->[vcf,tbi]}
			.collect(sort:true)
			.map{files->[[id:metadata.id],files]}
		)
	versions = versions.mix(BCFTOOLS_CONCAT.out.versions)

	VCF_TO_BED(
		BCFTOOLS_CONCAT.out.vcf.map{meta,vcf,tbi->[meta,vcf]}
		)
	versions = versions.mix(VCF_TO_BED.out.versions)

	BEDTOOLS_MERGE1(
		VCF_TO_BED.out.bed
		)
	versions = versions.mix(BEDTOOLS_MERGE1.out.versions)

	BED_CLUSTER(
		fasta,
		fai,
		dict,
		BEDTOOLS_MERGE1.out.bed
		)
	versions = versions.mix(BED_CLUSTER.out.versions)

	cluster_ch = BED_CLUSTER.out.bed
		.map{_meta,bed->bed}
		.map{it instanceof List?it:[it]}
		.flatMap()
		.map{bed->[[id:bed.baseName],bed]}
	
	BEDTOOLS_SLOP(
		fai,
		cluster_ch
		)
	versions = versions.mix(BEDTOOLS_SLOP.out.versions)

	BEDTOOLS_MERGE2(
		BEDTOOLS_SLOP.out.bed
		)
	versions = versions.mix(BEDTOOLS_MERGE2.out.versions)

emit:
    versions
	multiqc
	bed = BEDTOOLS_MERGE2.out.bed
	vcf = BCFTOOLS_CONCAT.out.vcf
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
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.charset.*;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.*;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;
import htsjdk.variant.vcf.VCFReader;

import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.*;


public class Minikit extends Main{
    private static final String TMP_FILE1 = "jeter1.vcf.gz";
    private static final String TMP_FILE2 = "jeter2.vcf.gz";
    private static final String TMP_FILE3 = "jeter3.vcf.gz";
    private static final Logger LOG = LogManager.getLogger(Main.class);
    //@Argument(fullName = StandardArgumentDefinitions.INTERVALS_LONG_NAME, shortName = StandardArgumentDefinitions.INTERVALS_SHORT_NAME, doc = "interval chr:start-end", common = true, optional = false)
    // private String interval = null;
    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "Variant input", common = true, optional = false)
	private String vcfIn;
    @Argument(fullName = StandardArgumentDefinitions.REFERENCE_LONG_NAME, shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME, doc = "Indexed fasta reference", common = true, optional = false)
	private String reference=null;
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "File to which variants should be written")
	private String finalVcf = null;
    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, doc = "A file containing the path to the bams", common = true, optional = false)
    private String bamsList;
    @Argument(fullName = "minimum-mapq", shortName = "Q", doc = "minimum mapping quality", common = true, optional = true)
	private int mapq= -1;
    @Argument(fullName = "report", shortName = "report", doc = "report", common = true, optional = false)
	private String reportOut="";

	private static class SimpleTimer {
		final long _start = System.currentTimeMillis();
		@Override
		public String toString() {
			final long d = System.currentTimeMillis() - this._start;
			return String.format("%.2f",d/1000.0);
			}
		}	


	private static class EqualRangeIteraror extends AbstractIterator<List<VariantContext>> {
		private final PeekableIterator<VariantContext> delegate;
		EqualRangeIteraror(CloseableIterator<VariantContext> delegate) {
			this.delegate= new PeekableIterator<>(delegate);
			}
		@Override
		protected List<VariantContext> advance() {
			if(!delegate.hasNext()) return null;
			final List<VariantContext> L = new ArrayList<>();
			final VariantContext ctx1 = delegate.next();
			L.add(ctx1);
			while(delegate.hasNext()) {
				final VariantContext ctx2 = delegate.peek();
				final int i = compare(ctx1,ctx2);
				if(i>0) throw new IllegalStateException(""+ctx1 + " " + ctx2);
				if(i<0) break;
				L.add(delegate.next());//consumme
				}
			return L;
			}
	
		}
	
	private static int compare(final VariantContext ctx1, final VariantContext ctx2) {
		int i= ctx1.getContig().compareTo(ctx2.getContig());
		if(i!=0) return i;
		return Integer.compare(ctx1.getStart(), ctx2.getStart());
		}
	
	private static String toString(final VariantContext ctx) {
		return ctx.getContig()+":"+ctx.getStart()+":"+ ctx.getAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining("/")) + "["
			+  ctx.getGenotypes().stream().filter(G->!(G.isNoCall() || G.isHomRef())).map(G->G.getSampleName()+":" + G.getGQ() + ":"+ G.getDP()+":"+G.getType().name()).collect(Collectors.joining(" "))
			+ "]";
		}

	private int diff(Path vcf1, Path vcf2,Path outvcfpath,int nExpect) throws IOException {
		final List<VariantContext> debug_list = new ArrayList<>();
		int count = 0;
		int discarded = 0;
		try(VCFReader vcfReader1  = new VCFFileReader(vcf1,false);
			VCFReader vcfReader2  = new VCFFileReader(vcf2,false);) {
			final VCFHeader header1 = vcfReader1.getHeader();
			final CloseableIterator<VariantContext> iter1 = vcfReader1.iterator();
			final CloseableIterator<VariantContext> iter2 = vcfReader2.iterator();
			final EqualRangeIteraror eq1 = new EqualRangeIteraror(iter1);
			final EqualRangeIteraror eq2 = new EqualRangeIteraror(iter2);

			
			
			final VariantContextWriterBuilder vcb = new VariantContextWriterBuilder();
			final SAMSequenceDictionary dict1 = header1.getSequenceDictionary();
			vcb.setReferenceDictionary(header1.getSequenceDictionary());
        	vcb.setIndexCreator(new TabixIndexCreator(dict1, TabixFormat.VCF));
        	vcb.setOption(Options.INDEX_ON_THE_FLY);
        	vcb.setOutputPath(outvcfpath);
        	vcb.setCreateMD5(false);
        	
        	final VariantContextWriter w = vcb. build();
        	w.writeHeader(header1);
        	
			while(eq1.hasNext()) {
				final List<VariantContext> array1 = new ArrayList<>(eq1.next());
				final VariantContext first1 = array1.get(0);
				
				List<VariantContext> array2 = null;
				while(eq2.hasNext()) {
					array2 = eq2.peek();
					final VariantContext first2 = array2.get(0);
					final int i= compare(first1, first2);
					if(i < 0) {
						array2 = null;
						break;
						}
					else if(i > 0) {
						eq2.next();//consumme && continue;
						continue;
						}
					break;
					}
				for(VariantContext ctx1 : array1) {
					boolean save=true;
					
					if(array2!=null) {
						final Set<Allele> alts = new HashSet<>(ctx1.getAlternateAlleles());
						alts.remove(Allele.SPAN_DEL);
						for(VariantContext ctx2 : array2) {
							if(compare(ctx1,ctx2)!=0) continue;
							if(!ctx1.getReference().equals(ctx2.getReference())) continue;
							ctx2.getGenotypes().stream().
								filter(G->!(G.isHomRef() || G.isNoCall())).
								flatMap(G->G.getAlleles().stream()).
								forEach(A->alts.remove(A));
							}
						if(alts.isEmpty()) save=false;
						}
					if(save) {
						w.add(ctx1);
						count++;
						if( debug_list.size()<10) debug_list.add(ctx1);
						}
					else
						{
						discarded++;
						LOG.info("remove "+toString(ctx1)+" was in " + array2.stream().map(V->toString(V)).collect(Collectors.joining(",")));
						}
					}
				}
			w.close();
			iter1.close();
			iter2.close();
		}
	if( count+ discarded != nExpect ) {
		throw new IllegalStateException("expected " + count + "+" + discarded + " == "+ nExpect);
		}
	LOG.info(getCommandLineName()+" : remains:  " + count + " discarded: "+ discarded);
	if(debug_list.size() < 10) LOG.info("Remain "+debug_list.stream().map(V->toString(V)).collect(Collectors.joining(" ")));
	return count;
	}

private static int count(final Path invcf) throws IOException {
	try(VCFReader vcfReader1  = new VCFFileReader(invcf,false)) {
		try(CloseableIterator<VariantContext> iter1 = vcfReader1.iterator()) {
			return (int)iter1.stream().count();
		}
	}
}
	
	
private static int copyTo(final Path invcf,final Path outvcfpath) throws IOException {
	final SimpleTimer timer = new SimpleTimer();
	int count=0;
	try(VCFReader vcfReader1  = new VCFFileReader(invcf,false)) {
		final VCFHeader header1 = vcfReader1.getHeader();
		
		
		final VariantContextWriterBuilder vcb = new VariantContextWriterBuilder();
			final SAMSequenceDictionary dict1 = header1.getSequenceDictionary();
			vcb.setReferenceDictionary(header1.getSequenceDictionary());
        	vcb.setIndexCreator(new TabixIndexCreator(dict1, TabixFormat.VCF));
        	vcb.setOption(Options.INDEX_ON_THE_FLY);
        	vcb.setOutputPath(outvcfpath);
        	vcb.setCreateMD5(false);
        	
        	try(final VariantContextWriter w = vcb. build()) {
	        	w.writeHeader(header1);
	        	try(final CloseableIterator<VariantContext> iter1 = vcfReader1.iterator()) {
	        		while(iter1.hasNext()) {
	        			w.add(iter1.next());
	        			count++;
	        			}
	        		}
        		}
		}
	LOG.info("Copy : "+invcf+" -> "+outvcfpath+" = "+count+" "+ timer+" seconds");
	return count;
	}



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


private int hc(Path tmpDir, final List<String> bams, final int nExpect, final BufferedWriter reportw, final SimpleTimer mainTimer) throws Exception {
	final SimpleTimer timer	= new SimpleTimer();
	LOG.info(getCommandLineName()+" " + bams.size()+ " bams : " + String.join(", ",bams));
	
	final Path in = tmpDir.resolve(TMP_FILE1);
	final Path in_tbi = tmpDir.resolve(TMP_FILE1 +".tbi");

	final Path out = tmpDir.resolve(TMP_FILE2);
	final Path out_tbi = tmpDir.resolve(TMP_FILE2 + ".tbi");
	
	final Path out2 = tmpDir.resolve(TMP_FILE3);
	final Path out2_tbi = tmpDir.resolve(TMP_FILE3 + ".tbi");

	
	final List<String> cmd= new ArrayList<>();
	cmd.add("HaplotypeCaller");
	cmd.add("-R");
	cmd.add(this.reference.toString());
	for(String bam: bams) {
		cmd.add("-I");
		cmd.add(bam);
		}
	cmd.add("-L");
	cmd.add(in.toString());
	cmd.add("--sample-ploidy");
	cmd.add("2");
	cmd.add("--do-not-run-physical-phasing");
	cmd.add("--alleles");
	cmd.add(in.toString());
	cmd.add("-O");
	cmd.add(out.toString());
	cmd.add("--tmp-dir");
	cmd.add(tmpDir.toString());
	if(this.mapq>0) {
		cmd.add("--minimum-mapping-quality");
		cmd.add(String.valueOf(this.mapq));
		}
	execute(cmd);
	
	int n = diff(in,out,out2, nExpect);
	
	Files.delete(in);
	Files.delete(in_tbi);
	Files.delete(out);
	Files.delete(out_tbi);
	Files.move(out2, in);
	Files.move(out2_tbi, in_tbi);

	LOG.info("That took "+ timer+" seconds.");
	reportw.write(String.valueOf(bams.size())+"\t"+String.join(",",bams) + "\t" + nExpect + "\t"+ n + "\t"+mainTimer+"\\n");
	return n;
}

private int doWork(final String[] args) {
	try {
		final  CommandLineArgumentParser cmdLineParser  = new CommandLineArgumentParser(this);
	        final boolean ret = cmdLineParser.parseArguments(System.err, args);
        	if (!ret) {
        		System.err.println(cmdLineParser.usage(false, false));
	        	return -1;
        		}

		final BufferedWriter w = Files.newBufferedWriter(Paths.get(this.reportOut),Charset.defaultCharset());
		//final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(Paths.get(this.reference));
		//final GenomeLocParser parser = new GenomeLocParser(dict);
		//final Locatable region  = parser.parseGenomeLoc(this.interval);
		
		final List<String> bams = new ArrayList<>(Files.readAllLines(Paths.get(this.bamsList)));

		if (bams.isEmpty()) {
	        	LOG.error("No Bam was provided");
        		return -1;
        		}
		final SimpleTimer timer= new SimpleTimer();
		Path tmpDir = Files.createTempDirectory("tmp");
		Path tmpVcf = tmpDir.resolve(TMP_FILE1);
		// copy original
		int nExpect = copyTo(Paths.get(vcfIn),tmpVcf);
		int nIter=0;
		while(!bams.isEmpty() && nExpect>0) {
			int n_bam_max = Math.min(50,nIter+1) ;
			if(nExpect>1000) n_bam_max=1;
			if(nExpect<100) n_bam_max=10;
			if(nExpect<10) n_bam_max=50;

			LOG.info(getCommandLineName() + " : remains "+bams.size()+" bams");
			final List<String> L = new ArrayList<>();
			L.add(bams.remove(0));
			while(!bams.isEmpty() && L.size() < n_bam_max) {
				L.add(bams.remove(0));
				}
			nExpect = hc(tmpDir,L , nExpect,w , timer);
			nIter++;
			}
		copyTo(tmpVcf,Paths.get(this.finalVcf));

		w.flush();
		w.close();
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

#####
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">compile gatk</entry>
</properties>
EOF
"""
}


process PER_VCF {
afterScript "rm -rf TMP"
tag "${bams} / ${vcf}"
memory "10g"
cpus 1
input:
	val(meta)
	val(genomeId)
	path(minikit)
	path(bams)
	path(vcf)
output:
	path("selection.vcf.gz"),emit:output
	path("report.tsv"),emit:report
	path("version.xml"),emit:version
script:
	def reference = params.genomes[genomeId].fasta
"""
hostname 1>&2
${moduleLoad("openjdk/11.0.8 bcftools")}

mkdir -p TMP
bcftools view -O z -o TMP/input.vcf.gz "${vcf.toRealPath()}" 
bcftools index -tf TMP/input.vcf.gz

java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -cp ${meta.gatkjar}:${minikit} Minikit \
                -I ${bams} \
                ${meta.mapq && ((meta.mapq as int)>0)?"--minimum-mapq ${meta.mapq}":""} \
                --output TMP/selection.vcf.gz \
                -V TMP/input.vcf.gz \
                --reference "${reference}" \
		--report TMP/report.tsv 


mv TMP/selection.vcf.gz ./
mv TMP/selection.vcf.gz.tbi ./
mv TMP/report.tsv ./

#####
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">find rare variants</entry>
</properties>
EOF
"""
}
