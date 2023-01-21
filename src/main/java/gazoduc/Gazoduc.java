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
package gazoduc;
import java.util.*;
import java.util.regex.*;
import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.function.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;


public class Gazoduc {
	public static final String ANSI_ESCAPE = "\u001B[";
	public static final String ANSI_RESET = ANSI_ESCAPE+"0m";
	private static final Logger LOG = Logger.getLogger(Gazoduc.class.getSimpleName());
	private static final String PARAM_GENOMES ="genomesFile";
	private static final String PARAM_GENOME ="genomeId";
	private static final String PARAM_REFERENCE="reference";
	private static final String MAIN_MENU="Main options";
	private static final String DEFAULT_VALUE= "value";
	public static final String DESC_INDEXED_FASTA = "Path to the reference genome as FASTA. The file must be indexed with 'samtools faidx' and 'samtools dict' ( or picard CreateSequenceDictionary ) ";
	public static final String DESC_INDEXED_BAM = "The BAM must be indexed with 'samtools index' (an associated .bai must be present) ";
	public static final String DESC_INDEXED_VCF = "The VCF must be indexed with 'bcftools index' (an associated .tbi/.csi must be present) ";
	public static final String DESC_VCF_OR_VCF_LIST = "Path to a VCF file or a file with the .list' suffix containing the full path to several VCFs file. " + DESC_INDEXED_VCF  ;
	public static final String DESC_VCF_LIST = "File with the .list' suffix, containing the full path to several VCFs file. " + DESC_INDEXED_VCF  ;
	public static final String DESC_BAM_OR_CRAM_LIST = "File with the .list' suffix, containing the full path to several BAMs ir CRAM file. " + DESC_INDEXED_BAM  ;
	public static final String DESC_JVARKIT_PEDIGREE = "Jvarkit formatted pedigree. Tab delimited, no header, FAM/ID/FATHER/MOTHER/SEX/PHENOTYPE";
	
	private static Gazoduc INSTANCE = null;
	private final Map<String,Object> params;
	private final List<Parameter> parameters = new Vector<>();
	private Genomes genomes = null;
	
	/**
	 * 
	 * wrapper for nextflow parameter
	 *
	 */
	public class Parameter {
		private final List<Validator> validators = new Vector<>();
		private final String key;
		private Object value = null;
		private String argName = DEFAULT_VALUE ;
		private String shortDesc = "";
		private String menu = MAIN_MENU;
		private boolean required = false;
		private boolean hidden = false;

		/**
		 * validator for parameter
		 * */
		private class Validator {
			public boolean validate() {
				final Object o = getParams().getOrDefault(getKey(),null);
				if(o==null) return true;
				return validateObject(o);
				}
			public boolean validateObject(final Object o) {
				if(o==null) return true;
				return validate(String.valueOf(o));
				}
			public boolean validate(final String o) {
				return true;
				}
			}
		/** constructor --key and value */
		private Parameter(final String key, Object value) {
			this.key = key;
			if(key==null || key.trim().isEmpty()) throw new IllegalArgumentException("Key is empty/null");
			if(key.contains("-"))  throw new IllegalArgumentException("Key contains '-' ("+key+")");
			this.value = value;
			}
		private Parameter(final String key) {
			this(key,null);
			}
		
		/** set value */
		public Parameter value(final Object o) {
			this.value = o;
			return this;
			}
		public Parameter argName(final String s) {
			this.argName = s;
			return this;
			}
		public Parameter menu(final String s) {
			this.menu = s;
			return this;
			}
		public Parameter desc(final String s) {
			this.shortDesc = s;
			return this;
			}

		public final Parameter description(final String s) {
			return this.desc(s);
			}

		public String getKey() {
			return this.key;
			}

		/* NF changes some keys when lower/uppercase change e.g: bwaPath -> bwa_path or selectAD -> select-AD */
		public boolean hasKey(final String s) {
			final String s1 = this.getKey().replaceAll("[\\-]","");
			final String s2 = s.replaceAll("[\\-]","");
			return s1.equalsIgnoreCase(s2);
			}

		public Parameter hidden() {
			this.hidden = true;
			return this;
			}
		public boolean isHidden() {
			return this.hidden;
			}

		public Parameter required() {
			this.required = true;
			return validator(new Validator() {
				@Override
				public boolean validate() {
					if(!getParams().containsKey(getKey()) || getParams().get(getKey())==null) {
						LOG.severe("option --"+ getKey() +" is not defined. " + markdown());
						return false;
						}
					return true;
					}
				});
			}

		public Parameter notEmpty() {
			return validator(new Validator() {
				@Override
				public boolean validate(final String s) {
					if(s==null || s.trim().isEmpty()) {
						LOG.severe("option --"+ getKey() +" is not defined or empty. " + markdown());
						return false;
						}
					return true;
					}
				});
			}

		
		private Parameter file(final Predicate<File> consummer, final String msg) {
			if(this.argName.equals(DEFAULT_VALUE)) argName("path to file");
			return notEmpty().
				validator(new Validator() {
				@Override
				public boolean validate(final String s) {
					try {
						final File f = new File(s);
						if(!consummer.test(f)) {
							LOG.severe("option --"+ getKey() +" = ("+s+") "+msg); 
							return false;
							}
						}
					catch(Throwable err) {
						LOG.severe("option --"+ getKey() +" = ("+s+") "+msg); 
						return false;
						}
					return true;
					}
				});
			}
		

		public Parameter existingFile() {
			if(this.argName.equals(DEFAULT_VALUE)) argName("path to file");
			return file(F->F.exists() && !F.isDirectory(),"Path should exist and be a file.");
			}

		public Parameter existingDirectory() {
			if(this.argName.equals(DEFAULT_VALUE) || this.argName.equals("file") || this.argName.equals("path")) {
				argName("path to directory");
				}
			return file(F->F.exists() && F.isDirectory(),"Path should exist and be a directory.");
			}
		
		public Parameter fullPath() {
			return file(F->F.exists() && F.toString().startsWith(File.separator),"File Should be full path");
			}
		
		/**
		 * check target is a reference is indexed with BWA
		 * @return this
		 */
		public Parameter bwaReference() {
			return desc("path to a reference indexed with bwa").
				file(F->{
					final File f = new File(F.getParentFile(), F.getName()+".bwt");
					return f.exists();
					},"Path should be the path to a fasta reference indexed with bwa");
			}

		public Parameter suffixes(final String...suffixes) {
			return validator(new Validator() {
				@Override
				public boolean validate(final String s) {					
					for(final String suffix : suffixes) {
						if(s.endsWith(suffix)) return true;	
						}
					LOG.severe("option --"+ getKey() +" value "+s+ "doesn't end with  (" + String.join(" | ", suffixes) +")/"); 
					return false;
					}
				});
			}

		public Parameter indexedFasta() {
			if(this.argName.equals(DEFAULT_VALUE)) argName("path to fasta file");
			return existingFile().
				suffixes(".fa",".fasta",".fa.gz",".fasta.gz",".fna",".fna.gz").
				file(F->{
						final String s= F.getName();
						File f = new File(F.getParentFile(),s+".fai");
						if(!f.exists())  {
							LOG.warning("not found "+f);
							return false;
							}
						int dot = s.lastIndexOf(".");
						if(dot<0) return false;
						String s2 = s.substring(0,dot);
						f = new File(F.getParentFile(),s2+".dict");
						if(!f.exists()) {
							LOG.warning("not found "+f);
							return false;
							}
						return true;
						}
					,"Fasta file must exists. ref.fa.fai and ref.dict must exists"
					);
				};
						
		public Parameter regex(final String reg) {
			return validator(new Validator() {
				@Override
				public boolean validate(final String s) {
					 Pattern p = Pattern.compile(reg);
					 Matcher m = p.matcher(s);
					 if(!m.matches()) {
						LOG.severe("option --"+ getKey() +" value "+s+ "doesn't match the regular expression '" + reg+"'"); 
						return false;
						}
					return true;
					}
				});
			}


		public Parameter inSet(final String...values) {
			return validator(new Validator() {
				@Override
				public boolean validate(final String s) {					
					for(final String s2 : values) {
						if(s.equals(s2)) return true;	
						}
					LOG.severe("option --"+ getKey() +" value "+s+ "doesn't end match with any of (" + String.join(" | ", values) +")/"); 
					return false;
					}
				});
			}

		public Parameter setBoolean() {
			return inSet("true","false").argName("true|false");
			}


		private Parameter setConsummer(java.util.function.Consumer<String> fun,final String msg) {
			return validator(new Validator() {
				@Override
				public boolean validate(final String s) {
					try {
						fun.accept(s);
						return true;
						}
					catch(final Throwable err) {
						LOG.severe("option --"+ getKey() +" value "+s+ " error. "+msg); 
						return false;
						}
					}
				});
			}

		/** add a validator to check the value is a {@link java.lang.Integer} */
		public Parameter setInteger() {
			if(this.argName.equals(DEFAULT_VALUE)) argName("integer");
			return setConsummer(S->Integer.parseInt(S),"Value should be an integer");
			}
		/** alias of setInteger */
		public final Parameter setInt() {
			return setInteger();
			}
		/** add a validator to check the value is a {@link java.lang.Long} */
		public Parameter setLong() {
			if(this.argName.equals(DEFAULT_VALUE)) argName("long");
			return setConsummer(S->Long.parseLong(S),"Value should be a long integer.");
			}
		/** add a validator to check the value is a {@link java.lang.Double} */
		public Parameter setDouble() {
			if(this.argName.equals(DEFAULT_VALUE)) argName("double");
			return setConsummer(S->Double.parseDouble(S),"Value should be a floating value.");
			}
		/** add a validator to check the value is a {@link URL} */
		public Parameter setURL() {
			if(this.argName.equals(DEFAULT_VALUE)) argName("url");
			return setConsummer(S->{
					try {
						new java.net.URL(S);
					} catch(MalformedURLException err) {
						throw new IllegalArgumentException(err);
					}
				},"Value should be an URL value.");
			}
		
		
		public Parameter setFile() {
			if(this.argName.equals(DEFAULT_VALUE)) argName("path");
			return setConsummer(S-> new File(S) ,"Value should be a file.");
			}

		public Parameter validator(final Validator v) {
			if(v!=null) this.validators.add(v);
			return this;
			}
		

		public boolean put() {
			if(getParams().containsKey(this.key)) { /** user set this key on the NF command line, it already defined but we register it */
				this.value = getParams().get(this.key);
				}
			else
				{
				getParams().put(this.key,this.value);
				}
			final Parameter old = Gazoduc.this.findParameterByName(this.key).orElse(null);
			if(old!=null) {
				//LOG.severe("key already defined in gazoduc ! "+ this.key+" " + old);
				return false;
				}
			else
				{
				Gazoduc.this.parameters.add(this);
				return true;
				}
			}



		public boolean validate() {
			boolean ok=true;
			for(Validator v: this.validators) {
				if(!v.validate()) ok=false;
				}
			if(!ok) {
				System.err.println("  --"+getKey()+"="+ this.value +" "+red("[FAILED]"));
				}
			else	{
				System.err.println("  --"+getKey()+"="+ this.value +" "+green("[OK]"));
				}
			return ok;
			}

		public String markdown() {
			return markdown(2);
			}
		public String markdown(int margin) {
			StringBuilder sb = new StringBuilder();
			while(margin>0) {
				sb.append(" ");
				margin--;
				}
			sb.append("* --");
			sb.append(this.key);
			sb.append(" <");
			sb.append(this.argName);
			sb.append("> . ");
			sb.append(this.shortDesc);
			sb.append(". ");
			if(this.required) {
				sb.append("[REQUIRED]. ");
				}
			if(this.value!=null) {
				sb.append("[");
				sb.append(this.value);
				sb.append("]");
				}
			sb.append(".\n");
			return sb.toString();
			}
		public String log() {
			return  new StringBuilder().append(this.key).append("\t|\t").append(this.value).toString();
			}
		@Override
		public String toString() {
			return markdown();
			}
		}

	public Parameter make(final String key, final Object value) {
		return new Parameter(key,value);
		}
	/** alias of make */
	public final Parameter build(final String key, final Object value) {
		return make(key,value);
		}
	/** 
	creates a new entry for a Parameter but don't insert it in the  <code>params</code> object. 
	 @param key name for this parameter
	 @param value value for this parameter
	 @return the new parameter
	 */
	public final Parameter param(final String key, final Object value) {
		return make(key,value);
		}

	/** 
	 alias of make with null value
	 @param key for this parameter
	 @return the new parameter
	 */
	public Parameter make(final String key) {
		return make(key,null);
		}
	/** alias of make
	 @return the new parameter
	 */
	public final Parameter build(final String key) {
		return make(key);
		}
	/** alias of make
	 @return the new parameter
	 */
	public final Parameter param(final String key) {
		return make(key);
		}
	
	/** find parameter by name
	 @return an {@link Optional} Parameter 
	 */
	public Optional<Parameter> findParameterByName(final String key) {
		return this.parameters.stream().filter(P->P.hasKey(key) ).findFirst();
		}

	public Gazoduc putGenomes() {
		make(PARAM_GENOMES,false).
			desc("Path to a XML file describing all the available genomes on your server. See doc").
			menu("Genomes").
			existingFile().
			required().
			put();
		
		make(PARAM_GENOME,false).
			desc( "The main genome used. This is the genome id in the XML file (see option --"+ PARAM_GENOMES+")").
			menu("Genomes").
			notEmpty().
			required().
			put();
		return this;
		}

	/**
	 * put the default parameters in the context. Default parameters are: <code>--help</code>, <code>--prefix</code> and  <code>--publishDir</code>
	 * @return this
	 */
	public Gazoduc putDefaults() {
		make("help",false).desc("Display help for this workflow and exit").menu("Help").setBoolean().put();
		make("prefix","").argName("string").desc("set a suffix for the files generated for this workflow").menu("Output").notEmpty().regex("[A-Z0-9a-z_\\.\\-]+").put();
		make("publishDir","").
			argName("directory").
			existingDirectory().
			fullPath().
			desc("set a base directory where final output files should be written.").
			menu("Output").
			notEmpty().
			put();
		return this;
		}
	
	/**
	 * put the default parameter  <code>>--reference</code> in the context
	 * @return this
	 */
	public Gazoduc putReference() {
		reference().put();
		return this;
		}

    public Gazoduc putGnomad() {
		make("gnomad_exome_hg19","/LAB-DATA/BiRD/resources/species/human/broadinstitute.org/gnomad/release-181127/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.vcf.gz").desc("gnomad exome path for hg19").menu("Gnomad").put();
        make("gnomad_genome_hg19","/LAB-DATA/BiRD/resources/species/human/broadinstitute.org/gnomad/release-181127/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.vcf.gz").desc("gnomad genome path for hg19").menu("Gnomad").put();
        make("gnomad_genome_hg38","/LAB-DATA/BiRD/resources/species/human/broadinstitute.org/gnomad/3.0/gnomad.genomes.r3.0.sites.vcf.gz").desc("gnomad genome path for hg38").menu("Gnomad").put();
		return this;
		}


    /**
   	 creates a parameter for <code>--reference</code> but doesn't put in the context
     @return the pararemeter for reference
     */
	public Parameter reference() {
		return make(PARAM_REFERENCE,false).
			argName("path to fasta").
			desc(DESC_INDEXED_FASTA).
			menu("Input").
			existingFile().
			required().
			indexedFasta();
		}

	public Gazoduc putCondaEnv() {
		final String env = System.getenv("CONDA_ENVS_PATH");
		make("conda",env).
			argName("CONDA_ENVS_PATH").
			desc("The bird Cluster at institut du Thorax as a global env variable ${CONDA_ENVS_PATH} where conda finds its file").
			menu("Conda").
			notEmpty().
			required().
			put();
		return this;
		}
		
	public class UsageBuilder {
		private String name = "workflow";
		private String description = "no description available";
		private List<String> authors = new ArrayList<>();
		private UsageBuilder() {
			this.authors.add("Pierre Lindenbaum PhD. Institut du Thorax. U1087. 44400 Nantes. France.");
			}

		public UsageBuilder name(final String s) {
			this.name = s;
			return this;
			}
		public UsageBuilder description(final String s) {
			this.description = s;
			return this;
			}

		public UsageBuilder desc(final String s) {
			return description(s);
			}

		public UsageBuilder author(final String s) {
			authors.add(s);
			return this;
			}
		
		public void print() {
			System.out.println(markdown());
			}

		public void log() {
			final StringBuilder sb = new StringBuilder();
			for(Parameter p: Gazoduc.this.parameters) {
				if(p.isHidden()) continue;
				sb.append(p.log()).append("\n");
				}
			LOG.info(sb.toString());
			}

		
		public String markdown() {
			final Set<String> menus = Gazoduc.this.parameters.stream().
				filter(S->!S.isHidden()).
				map(S->S.menu).
				collect(Collectors.toCollection(TreeSet::new));

			final StringBuilder w = new StringBuilder();

			w.append("# "+ this.name+"\n\n");
			w.append(this.description+"\n\n");

			w.append("## Author(s)\n\n");
			
			for(String a: this.authors) {
				w.append("  * ").append(a).append("\n");
				}
			w.append("\n");
			
			w.append("## Options\n\n");
			
			for(int side=0;side<2;side++) {
			for(String menu: menus) {
				if(menu.equals(MAIN_MENU) && side==1) continue;
				if(!menu.equals(MAIN_MENU) && side==0) continue;
				w.append("### ").append(menu.isEmpty()?"Main options":menu).append("\n\n");
				final String fmenu = menu;
				for(Parameter p: Gazoduc.this.parameters.stream().
							filter(S->S.menu.equals(fmenu)).
							sorted((A,B)->A.key.compareTo(B.key)).
							collect(Collectors.toList())) {
					w.append(p.markdown(4));
					}
				w.append("\n");
				}
			w.append("\n");
			}
			w.append("## Issues\n\n");
			w.append("report issues at https://github.com/lindenb/gazoduc-nf/issues\n\n");
			
			try {
				final File f = new File("workflow.svg");
				if(f.exists()) {
					w.append("## Workflow\n\n");
					w.append("![workflow](./workflow.svg)\n\n");
					}
				}
			catch(final Throwable err) {
				}
			
			return w.toString();
			}
		@Override
		public String toString() {
			return markdown();
			}
		}


	public UsageBuilder usage() {
		return new UsageBuilder();
		}
	

	public void validate() {
		System.err.println("VALIDATION");
		System.err.println("==========");
		boolean is_valid = true;
		for(Parameter p: this.parameters) {
			if(!p.validate()) is_valid=false;
			}
		for(String key : getParams().keySet()) {
			if( this.findParameterByName(key).isPresent() ) continue;
			System.err.println("key \"--"+key+"\" was defined in params but was not declared ["+yellow("WARNING")+"].");
			}
		if(getParams().containsKey(PARAM_GENOMES)) {
			try {
				final Genomes genomes  = this.getGenomes();
				 if(getParams().containsKey(PARAM_GENOME)) {
					try {
						if(!this.getGenome().validate()) {
							is_valid = false;
							}
						}
					 catch(Throwable err2) {
		                                System.err.println("Cannot load genome ["+red("ERROR")+"]." + err2.getMessage());
                		                is_valid = false;
                                		}
					}
				}
			catch(Throwable err) {
				System.err.println("Cannot load XML file for genomes ["+red("ERROR")+"]." + err.getMessage());
				is_valid = false;
				}
			}
		else if(getParams().containsKey(PARAM_GENOME)) {
			System.err.println("--"+PARAM_GENOME+" is defined but not --"+PARAM_GENOMES+" ["+yellow("WARNING")+"].");
			}


		if(!is_valid) {
			throw new IllegalArgumentException("Validation of parameters failed");
			}
		}

	public Genomes getGenomes() {
		if(this.genomes == null) {
			if(!findParameterByName(PARAM_GENOMES).isPresent()) {
				throw new IllegalArgumentException("Asking for genomes but --"+PARAM_GENOMES+" was not defined");
				}
			final String filename = getParams().getOrDefault(PARAM_GENOMES, "NO_FILE").toString();
			try {
				this.genomes = Genomes.load(filename);
			} catch (final Exception e) {
				LOG.log(Level.SEVERE,"Cannot load genomes from "+ filename, e);
				throw new RuntimeException("Cannot load genomes --"+PARAM_GENOMES+"="+filename, e);
				}
			}
		return this.genomes;
		}
	
	public Genome getGenome() {
		final Genomes genomes = getGenomes();
		Optional<Parameter> opt= findParameterByName(PARAM_GENOME);
		if(opt.isPresent()) {
			return genomes.getById(getParams().getOrDefault(PARAM_GENOME,"").toString());
			}
		 opt= findParameterByName(PARAM_REFERENCE);
		 if(opt.isPresent()) {
			return genomes.getByReference(getParams().getOrDefault(PARAM_REFERENCE,"").toString());
			}
		throw new RuntimeException("Cannot load current genome. No parameter  --"+PARAM_GENOME +" or --"+PARAM_REFERENCE+" is defined.");
		}
	
	private static String pen(int pen,String s) {
		return ANSI_ESCAPE+pen+"m"+s+ANSI_RESET;
		}

	private static String green(final String s) {
		return pen(32,s);
		}

	private static String red(final String s) {
		return pen(31,s);
		}

	private static String yellow(final String s) {
		return pen(93,s);
		}
	
	private Gazoduc(final Map<String,Object> params) {
		this.params = params;
		if(params==null) throw new IllegalArgumentException("params is null");
		}

	private Map<String,Object> getParams() {
		return this.params;
		}

	public static Gazoduc getInstance(final Map<String,Object> params) {
		if(INSTANCE==null) {
			synchronized(Gazoduc.class) {
				INSTANCE = new Gazoduc(params);
				}
			}
		return INSTANCE;
		}
	/** return current singleton instance */
	public static Gazoduc getInstance() {
		if(INSTANCE==null) throw new IllegalStateException("gazoduc.Gazoduc was not initialized. INSTANCE is null.");
		return INSTANCE;
		}
	}
