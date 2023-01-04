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
import java.util.function.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;


public class Gazoduc {
	private static final Logger LOG = Logger.getLogger(Gazoduc.class.getSimpleName());
	private static final String MAIN_MENU="Main options";
	private static final String DEFAULT_VALUE= "value";
	private static Gazoduc INSTANCE = null;
	private final List<Parameter> parameters = new Vector<>();



	public class Parameter {
		private final List<Validator> validators = new Vector<>();
		private final String key;
		private Object value = null;
		private String argName = DEFAULT_VALUE ;
		private String shortDesc = "";
		private String longDesc = "";
		private String menu = MAIN_MENU;
		private boolean required = false;

		private class Validator {
			public boolean validate(final Map<String,Object> map) {
				Object o = map.getOrDefault(getKey(),null);
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

		Parameter(final String key, Object value) {
			this.key = key;
			if(key==null || key.trim().isEmpty()) throw new IllegalArgumentException("Key is empty/null");
			this.value = value;
			}
		Parameter(final String key) {
			this(key,null);
			}
		public Parameter value(Object o) {
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
		public String getKey() {
			return this.key;
			}
		public Parameter required() {
			this.required = true;
			return validator(new Validator() {
				@Override
				public boolean validate(final Map<String,Object> map) {
					if(!map.containsKey(getKey()) || map.get(getKey())==null) {
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


		public Parameter existingFile() {
			if(this.argName.equals(DEFAULT_VALUE)) argName("path to file");
			return notEmpty().
				validator(new Validator() {
				@Override
				public boolean validate(final String s) {
					try {
						final File f = new File(s);
						if(!f.exists()) {
							LOG.severe("option --"+ getKey() +" file doesn't exist: "+f);
							return false; 
							}
						}
					catch(Throwable err) {
						LOG.severe("option --"+ getKey() +" is not defined or empty. " + markdown()+" "+err.getMessage()); 
						return false;
						}
					return true;
					}
				});
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
				validator(new Validator() {
                                @Override
                                public boolean validate(final String s) {
					boolean ok=true;
					try {
						File f = new File(s);
						if(!f.exists()) ok=false;
						f = new File(s+".fai");
						if(!f.exists()) ok=false;
						int dot = s.lastIndexOf(".");
						if(dot>0 ) {
							String s2 = s.substring(0,dot);
							f = new File(s2+".dict");
							if(!f.exists()) ok=false;
							}
						return true;
						}
					catch(Throwable err) {
						LOG.severe(err.getMessage());
						}
					if(!ok) LOG.severe("option --"+ getKey() +" value "+s+ "doesn't look like an indexed fasta sequence (.fai and .dict required)/");
                                        return ok;
                                        }
				});
			}

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
			return inSet("true","false");
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


		public Parameter setInteger() {
			if(this.argName.equals(DEFAULT_VALUE)) argName("integer");
			return setConsummer(S->Integer.parseInt(S),"Value should be an integer");
			}

		public Parameter setLong() {
			if(this.argName.equals(DEFAULT_VALUE)) argName("long");
			return setConsummer(S->Long.parseLong(S),"Value should be a long integer.");
			}

		public Parameter setDouble() {
			if(this.argName.equals(DEFAULT_VALUE)) argName("double");
			return setConsummer(S->Double.parseDouble(S),"Value should be a floating value.");
			}



		public Parameter validator(final Validator v) {
			if(v!=null) this.validators.add(v);
			return this;
			}
		

		public boolean put(final Map<String,Object> map) {
			if(map.containsKey(this.key)) { /** user set this key on the NF command line, it already defined but we register it */
				this.value = map.get(this.key);
				}
			else
				{
				map.put(this.key,this.value);
				}
			final Parameter old = Gazoduc.this.parameters.stream().filter(P->P.key.equals(this.key)).findFirst().orElse(null);
			if(old!=null) {
				LOG.severe("key already defined in gazoduc ! "+ this.key+" " + old);
				return false;
				}
			else
				{
				Gazoduc.this.parameters.add(this);
				return true;
				}
			}

		public boolean validate(final Map<String,Object> map) {
			boolean ok=true;
			for(Validator v: this.validators) {
				if(!v.validate(map)) ok=false;
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
				sb.append("[ ");
				sb.append(this.value);
				sb.append("]");
				}
			sb.append(".\n");
			return sb.toString();
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
	/** alias of make */
	public final Parameter param(final String key, final Object value) {
		return make(key,value);
		}


	public Parameter make(final String key) {
		return make(key,null);
		}
	/** alias of make */
	public final Parameter build(final String key) {
		return make(key);
		}
	/** alias of make */
	public final Parameter param(final String key) {
		return make(key);
		}
	


	public Gazoduc putDefaults(final Map<String,Object> map) {
		make("help",false).desc("Display help for this workflow and exit").menu("Help").setBoolean().put(map);
		make("prefix","").argName("string").desc("set a suffix for the files generated for this workflow").menu("Output").notEmpty().regex("[A-Z0-9a-z_\\.\\-]+").put(map);
		make("publishDir","").argName("directory").desc("set a base directory where final output files should be written.").menu("Output").notEmpty().put(map);
		return this;
		}
	public Parameter reference() {
		return make("reference",false).
			argName("path to fasta").
			desc("Path to the reference genome as FASTA. The file must be indexed with 'samtools faidx' and 'samtools dict'").menu("Input").
			indexedFasta();
		}

	
	public class UsageBuilder {
		private UsageBuilder() {
			}
		public void print() {
			System.out.println(markdown());
			}
		public String markdown() {
			final Set<String> menus = Gazoduc.this.parameters.stream().
				map(S->S.menu).
				collect(Collectors.toCollection(TreeSet::new));

			StringBuilder w = new StringBuilder();
			
			for(int side=0;side<2;side++) {
			for(String menu: menus) {
				if(menu.equals(MAIN_MENU) && side==1) continue;
				if(!menu.equals(MAIN_MENU) && side==0) continue;
				w.append("## ").append(menu.isEmpty()?"Main options":menu).append("\n\n");
				for(Parameter p: Gazoduc.this.parameters) {
					if(!p.menu.equals(menu)) continue;
					w.append(p.markdown(4));
					}
				w.append("\n");
				}
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
	

	public void validate(final Map<String,Object> m) {
		boolean is_valid = true;
		for(Parameter p: this.parameters) {
			if(!p.validate(m)) is_valid=false;
			}
		if(!is_valid) {
			throw new IllegalArgumentException("Validation failed");
			}
		}


	public static Gazoduc getInstance() {
		if(INSTANCE==null) {
			synchronized(Gazoduc.class) {
				INSTANCE = new Gazoduc();
				}
			}
		return INSTANCE;
		}

	
	}
