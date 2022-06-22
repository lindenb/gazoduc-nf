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

include {getModules} from '../utils/functions.nf'

process COMPILE_GTF_TO_BED {
executor "local"
afterScript "rm -rf TMP"
input:
	val(tag)
output:
	path("gtf2bed.jar"),emit:jar
	path("version.xml"),emit:version

script:
"""
module load ${getModules("java")}

mkdir TMP

cat << __EOF__ > TMP/GTFToBed.java
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.PushbackInputStream;
import java.io.PrintStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.HashSet;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

public class GTFToBed {
private	final Set<String> columns = new TreeSet<>();
private final Pattern tab = Pattern.compile("[\\t]");
private enum Type{ undefined,gtf,gff3}
private Type gtype = Type.undefined;

private static InputStream mayGZIP(InputStream inputStream) throws IOException {
        final PushbackInputStream pushbackInputStream = new PushbackInputStream(inputStream, 2);
        final byte[] signature = new byte[2];
        final int length = pushbackInputStream.read(signature);
        pushbackInputStream.unread(signature, 0, length);

        boolean isGzipped = ((signature[0] == (byte) (GZIPInputStream.GZIP_MAGIC)) && (signature[1] == (byte) (GZIPInputStream.GZIP_MAGIC >> 8)));
        if (isGzipped) {
            final int kb64 = 65536;
            return new GZIPInputStream(pushbackInputStream, kb64);
        } else {
            return pushbackInputStream;
        }
    }


private void put(Map<String,String> map,String key,String value) {
	if(key.isEmpty()) throw new IllegalArgumentException("key is empty ");
	if(!this.columns.contains(key)) return;
	if(map.put(key, value)!=null) {
		throw new IllegalArgumentException("duplicate key "+key+" in "+map);
		}
	}

private Map<String,String> attributes(final String s) {
	final Map<String,String> map= new HashMap<>();
	if(s.equals(".")) return map;
	int i=0;
	
	while(i< s.length()) {
		/* skip ws */
		while(i< s.length() && Character.isWhitespace(s.charAt(i))) i++;
		final StringBuilder key=new StringBuilder();
		while(i< s.length()) {
			if(this.gtype.equals(Type.undefined)) {
				if(Character.isWhitespace(s.charAt(i))) {
					this.gtype = Type.gtf;
					i++;
					break;
					}
				else if(s.charAt(i)=='=') {
					this.gtype = Type.gff3;
					i++;
					break;
					}
				}
			if(this.gtype.equals(Type.gtf) && Character.isWhitespace(s.charAt(i))) {
				i++;
				break;
				}
			else if(this.gtype.equals(Type.gff3) && s.charAt(i)=='=') {
				i++;
				break;
				}
			key.append(s.charAt(i));
			i++;
			}
		if(this.gtype.equals(Type.undefined)) throw new IllegalArgumentException("undefined gtf type with " + s);
		/* skip ws */
		while(i< s.length() && Character.isWhitespace(s.charAt(i))) i++;
		if(i>=s.length()) throw new IllegalArgumentException("expected '\\"' or '=' after "+s.substring(0,i));
		
		if(this.gtype.equals(Type.gtf)) {
			if(s.charAt(i)!='\\"')  throw new IllegalArgumentException("expected quote after "+s.substring(0,i));
			i++;
			}
		
		
		final StringBuilder value=new StringBuilder();
		while(i< s.length() ) {
			if(this.gtype.equals(Type.gtf) && s.charAt(i)=='\\"') break;
			else if(this.gtype.equals(Type.gff3) && s.charAt(i)==';') break;
			
			if(s.charAt(i)=='\\\\') {
				if(i+1>=s.length()) throw new IllegalArgumentException("unknown escape sequence after "+s.substring(0,i));
				i++;
				switch(s.charAt(i)) {
					case '"': value.append("\\\\"");break;
					case '\\'': value.append("\\'");break;
					case '\\\\': value.append("\\\\");break;
					case 't': value.append("\\t");break;
					case 'n': value.append("\\n");break;
					default: throw new IllegalArgumentException("unknown escape sequence after "+s.substring(0,i));
					}
				}
			else
				{
				value.append(s.charAt(i));
				}
			i++;
			}
		if(this.gtype.equals(Type.gtf)) {
			if(i>=s.length()) throw new IllegalArgumentException("expected quote after "+s.substring(0,i));
			if(this.gtype.equals(Type.gtf) &&  s.charAt(i)!='\\"')  throw new IllegalArgumentException("expected quote after "+s.substring(0,i));
			i++;
			}
		this.put(map,key.toString(), value.toString());
		
		/* skip ws */
		while(i< s.length() && Character.isWhitespace(s.charAt(i))) i++;
		if(i<s.length()) {
			if(s.charAt(i)!=';')  throw new IllegalArgumentException("expected semicolon after "+s.substring(0,i));
			i++;
			}
		/* skip ws */
		while(i< s.length() && Character.isWhitespace(s.charAt(i))) i++;
		}
	return map;
}

private void instanceMain(final String args[]) {
	try {
		final String MY_COLS="gtf.source,gtf.feature,gtf.score,gtf.strand,gtf.frame";

		final String GTF_COLS = "ccds_id,exon_id,exon_number,exon_version,gene_biotype,gene_id,gene_name,gene_source,"
			+ "gene_version,havana_transcript,havana_transcript_version,protein_id,protein_version,tag,"
			+ "transcript_biotype,transcript_id,transcript_name,transcript_source,transcript_version";
		final String GFF_COLS= 
			"Alias,biotype,ccdsid,constitutive,description,ensembl_end_phase,"
				+ "ensembl_phase,exon_id,external_name,gene_id,havana_transcript,"
				+ "havana_version,ID,logic_name,Name,Parent,protein_id,rank,tag,"
				+ "transcript_id,version"
				;
		String colStr = MY_COLS+","+GTF_COLS+","+GFF_COLS;
		int optind=0;
		while(optind < args.length) {
			if(args[optind].equals("--columns") && optind+1 < args.length) {
				optind++;
				colStr=args[optind];
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
		
		this.columns.clear();
		this.columns.addAll(Arrays.asList(colStr.split("[, \\t;\\\\|]+")));
		this.columns.remove("");
		
		PrintStream out = System.out;
		out.println("#chrom\\tstart\\tend"+String.join("\\t",this.columns));
		
		try(BufferedReader br = new BufferedReader(new InputStreamReader(mayGZIP(System.in), "UTF-8"))) {
			br.lines().
				filter(S->!S.startsWith("#")).
				map(S->tab.split(S)).
				forEach(tokens->{
					if(tokens.length!=9) throw new IllegalArgumentException("expected 9 tokens in "+String.join("<\\\\t>", tokens));
					out.print(tokens[0]);
					out.print('\\t');
					out.print(Integer.parseInt(tokens[3])-1);
					out.print('\\t');
					out.print(Integer.parseInt(tokens[4]));
					
					final Map<String,String> atts = attributes(tokens[8]);

					
					put(atts,"gtf.source",tokens[1]);
					put(atts,"gtf.feature",tokens[2]);
					put(atts,"gtf.score",tokens[5]);
					put(atts,"gtf.strand",tokens[6]);
					put(atts,"gtf.frame",tokens[7]);
					
					
					for(final String field:this.columns) {
						out.print('\\t');
						out.print(atts.getOrDefault(field,"."));
						}
					out.println();
				});
			out.flush();
			}
		}
	catch(final Throwable err ) {
		err.printStackTrace();
		System.exit(-1);
	}
}
public static void main(final String[] args) {
	new GTFToBed().instanceMain(args);
	}
}
__EOF__


cat < EOF > TMP/tmp.mf
Manifest-Version: 1.0
Main-Class: GTFToBed
EOF


javac -d TMP -sourcepath TMP TMP/GTFToBed.java
jar cfm gtf2bed.jar TMP/tmp.mf -C TMP .

###############################################################################
cat <<EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">compile gtf2bed</entry>
	<entry key="javac.version">\$(javac -version 2>&1)</entry>
</properties>
EOF
"""
}
