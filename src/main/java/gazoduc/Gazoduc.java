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
import java.util.function.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;


public class Gazoduc {
	private static final Logger LOG = Logger.getLogger(Gazoduc.class.getSimpleName());
	private static Gazoduc INSTANCE = null;
	private final List<Parameter> parameters = new Vector<>();


	public class Parameter {
		private final String key;
		private Object value;
		private String argName = "value";
		private String shortDesc = "";
		private String longDesc = "";
		private String menu = "";
		private boolean required = false;
		Parameter(final String key, Object value) {
			this.key = key;
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
		public Parameter required(final boolean b) {
			this.required = b;
			return this;
			}
		public Parameter required() {
			return this.required(true);
			}
		public void put(final Map<String,Object> map) {
			if(map.containsKey(this.key)) {
				LOG.error("params already contains "+this);
				return -1;
				}
			Gazoduc.this.parameters.add(this);
			map.put(this.key,this.value);
			}
		@Override
		public String toString() {
			}
		}

	public Parameter make(final String key, final Object value) {
		return new Parameter(key,value);
		}
	
		public String markdown() {
			StringBuilder sb = new StringBuilder(" * ");
			sb.append(getName());
			sb.append(" <");
			sb.append(getArgName());
			sb.append("> ");
			sb.append(getDescription());
			sb.append("\n");
			return sb.toString();
			}
		}
	
	/** parameter factory, creats a parameter */
	public class ParameterFactory extends BaseParam {
		ParameterFactory(final String name) {
			super(name);
			}
		ParameterFactory description(final String s) {
			return description(()->s);
			} 
		ParameterFactory description(final Supplier<String> s) {
			this.description = s;
			return this;
			} 
		ParameterFactory menu(final String s) {
			this.submenu = s;
			return this;
			}
		ParameterFactory argName(final String s) {
			this.argName = s;
			return this;
			}
		ParameterFactory validator(Validator v) {
			if(v!=null) validators.add(v);
			return this;
			}
		
		
		ParameterFactory setInteger() { return setClass(Integer.class);}
		ParameterFactory setLong() { return setClass(Long.class);}
		ParameterFactory setClass(final Class<?> C) {
			return validator(new Validator() {
				@Override
				public boolean validate(Parameter p,Map<String,Object> o) {
					return true;
					}
				});
			}
		ParameterFactory value(final String s) {
			return value(()->s);
			} 
		ParameterFactory value(final Supplier<String> s) {
			this.value = s;
			return this;
			} 
		public Parameter build(Map<String,Object> nfParams) {
			Parameter p = new Parameter(this.name);
			p.copy(this);
			Gazoduc.this.params.put(p.getName(),p);
			String v = p.value.get();
			if(v!=null) nfParams.put(p.getName(),v);
			return p;
			}
		}

	private Map<String,Parameter> params = new HashMap<>();

	private Gazoduc() {
		}
	
	public ParameterFactory option(final String name) {
		return new ParameterFactory(name);
		}

	public void validate(final Map<String,Object> m) {
		boolean is_valid = true;
		for(Parameter p: this.params.values()) {
			if(!p.validate(m)) is_valid=false;
			}
		if(!is_valid) {
			
			}
		}

	public String usage() {
		final Set<String> menus = this.params.values().stream().
				map(S->S.submenu).
				collect(Collectors.toCollection(TreeSet::new));

		StringBuilder w = new StringBuilder();

		for(String menu: menus) {
			w.append("## ").append(menu.isEmpty()?"Main options":menu).append("\n\n");
			for(Parameter p: this.params.values()) {
				if(!p.submenu.equals(menu)) continue;
				w.append(p.markdown());
				}
			}
		return w.toString();
		}

	public static Gazoduc getInstance() {
		if(INSTANCE==null) {
			synchronized(Gazoduc.class) {
				INSTANCE = new Gazoduc();
				}
			}
		return INSTANCE;
		}
	public String what(Object o) {	
		Class<?> c = o.getClass();
		do {
		    System.err.println(c);
			c=c.getSuperclass();
                    } while(!c.equals(Object.class));
		
		return String.valueOf(o);
		}

	public boolean parseBoolean(final Object o) {	
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
