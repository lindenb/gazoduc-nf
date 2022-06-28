package gazoduc;
import java.util.*;
import java.util.function.*;
import java.util.stream.Collectors;


public class Gazoduc {
	private static Gazoduc INSTANCE = null;

	public abstract class Validator {
		abstract boolean validate(Parameter p,Map<String,Object> o);
		}

	public abstract class BaseParam {
		final String name;
		String submenu = "";
		boolean help = false;
		Supplier<String> description = ()->this.getName();
		Supplier<String> value = ()->null;
		List<Validator> validators = new ArrayList<>();
		BaseParam(final String name) {
			this.name = name;
			}
		String getName() {
			return this.name;
			}
		BaseParam copy(BaseParam o) {
			if(o==this) return this;
			this.submenu = o.submenu;
			this.help = o.help;
			this.description = o.description;
			this.value = o.value;
			this.validators = new ArrayList<>(o.validators);
			return this;
			}
		}

	public class Parameter extends BaseParam {
		Parameter(final String name) {
			super(name);
			}
		String getArgName() {
			return "arg";
			}
		boolean validate(Map<String,Object> map) {
			boolean is_valid = true;
			for(Validator v: this.validators) {
				if(!v.validate(this,map)) is_valid =false;
				}
			return is_valid;
			}
		public String getDescription() {
			return this.description.get();
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
		ParameterFactory validator(Validator v) {
			if(v!=null) validators.add(v);
			return this;
			}
		ParameterFactory setInteger() { return setClass(Integer.class);}
		ParameterFactory setLong() { return setClass(Long.class);}
		ParameterFactory setClass(final Class<?> C) {
			return validator(new Validator() {
				@Override boolean validate(Parameter p,Map<String,Object> o) {
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

	public boolean parseBoolean(Object o) {	
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
