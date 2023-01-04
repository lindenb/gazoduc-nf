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
import java.io.*;
import java.util.function.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;


public class Gazoduc {
	private static final Logger LOG = Logger.getLogger(Gazoduc.class.getSimpleName());
	private static final String MAIN_MENU="Main options";
	private static Gazoduc INSTANCE = null;
	private final List<Parameter> parameters = new Vector<>();



	public class Parameter {
		private final List<Validator> validators = new Vector<>();
		private final String key;
		private Object value = null;
		private String argName = "value";
		private String shortDesc = "";
		private String longDesc = "";
		private String menu = MAIN_MENU;


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





		public Parameter validator(final Validator v) {
			if(v!=null) this.validators.add(v);
			return this;
			}
		

		public boolean put(final Map<String,Object> map) {
			if(map.containsKey(this.key)) {
				LOG.warning("params already contains "+this);
				return false;
				}
			Gazoduc.this.parameters.add(this);
			map.put(this.key,this.value);
			return true;
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
			sb.append("> ");
			sb.append(this.shortDesc);
			sb.append("\n");
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

	static boolean parseBoolean(final Object o) {	
		if(o==null) return false;
		String s= String.valueOf(o).toLowerCase();
		if(s.equals("t")) return true;
		if(s.equals("true")) return true;
		if(s.equals("yes")) return true;
		if(s.equals("y")) return true;
		if(s.equals("on")) return true;
		if(s.equals("1")) return true;


		if(s.equals("f")) return false;
		if(s.equals("false")) return false;
		if(s.equals("no")) return false;
		if(s.equals("n")) return false;
		if(s.equals("off")) return false;
		if(s.equals("0")) return false;

		throw new IllegalArgumentException("Illegal boolean value "+s+".");
		}
	
	}
