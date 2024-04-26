import java.io.BufferedReader;
import java.util.regex.*;
import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.stream.*;

public class Minikit {

private static class MinMax {
		Double min=null;
		Double max= null;
		}

private final Map<String,MinMax> tags = new HashMap<>();

private void fill(String key,String value) {
if(value.equals(".") || value.equalsIgnoreCase("nan")) return;
Double n;
try {
	n = new Double(value);
} catch(Throwable err) {
	return;
}
MinMax mm = tags.get(key);
if(mm.min==null) {
	mm.min= n;
	mm.max= n;
	}
else	{
	if(n.compareTo(mm.min)<0) mm.min=n;
	if(mm.max.compareTo(n)<0) mm.max=n;
	}
}

private void instanceMain(final String args[]) {
	try {
		final Pattern comma = Pattern.compile("[,]");
		final Pattern tab = Pattern.compile("[\t]");
		final Pattern semicolon = Pattern.compile("[;]");
		String mode="";
		int optind=0;
		while(optind < args.length) {
			if(args[optind].equals("--an") && optind+1< args.length) {
				optind++;
				for(String tag : args[optind].split("[, ;]+")) {
					tags.put(tag,new MinMax());
					}
				}
			else if(args[optind].equals("--")) {
				optind++;
				break;
				}
			else if(args[optind].startsWith("-")) {
				System.err.println("unknown option "+args[optind]);
				System.exit(-1);
				}
			optind++;
			}
		if(optind != args.length) {
			System.err.println("illegal number of arguments");
			System.exit(-1);
			}
		try(BufferedReader br = new BufferedReader(new InputStreamReader(System.in))) {
			br.lines().
				filter(L->!L.startsWith("#")).
				map(T->tab.split(T)).
				map(T->T[7]).
				flatMap(T->Arrays.stream(semicolon.split(T))).
				forEach(T->{
					final int i= T.indexOf("=");
					if(i==-1) return;
					final String key = T.substring(0,i);
					if(!tags.containsKey(key)) return;
					final String value = T.substring(i+1);
					if(value.contains(",")) return;
					fill(key,value);
					});
			}
		catch(Exception err) {
			throw err;
			}
		int n=0;
		PrintStream out = System.out;
		for(String key: tags.keySet()) {
			final MinMax mm = tags.get(key);
			if(mm.min==null) continue;
			if(mm.max==null) continue;
			if(mm.min.equals(mm.max)) continue;
			out.println("-an");
			out.println(key);
			n++;
			}
		out.flush();
		if(n==0) System.err.println("NO VALID TAGS WAS FOUND IN "+String.join(" , ",tags.keySet()));
	    }
	catch(final Throwable err ) {
		err.printStackTrace();
		System.exit(-1);
		}
	}
public static void main(final String[] args) {
	new Minikit().instanceMain(args);
	}
}
