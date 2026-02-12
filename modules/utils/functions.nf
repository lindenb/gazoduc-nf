/*

Copyright (c) 2026 Pierre Lindenbaum

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

import static org.fusesource.jansi.Ansi.Color 

//static java.util.Set<String> __warnings = new java.util.HashSet<>();


/** return boolean from hasmap */
def getBoolean(meta,key) {
	if(!meta.containsKey(key)) return false;
	return parseBoolean(meta.get(key));
	}


def parseBoolean(def o) {
	if(o == null)  {
		def msg =  "cannot convert null to boolean";
	        log.warn(msg);
        	throw new IllegalArgumentException(msg);
		}
	if((o instanceof Boolean)) {
		return (Boolean)o;
		}
	final String s= String.valueOf(o).trim().toLowerCase();
	if(s.equals("true")) return true;
	if(s.equals("yes")) return true;
	if(s.equals("on")) return true;
	if(s.equals("t")) return true;
	if(s.equals("y")) return true;
	if(s.equals("1")) return true;

	if(s.equals("false")) return false;
	if(s.equals("no")) return false;
	if(s.equals("off")) return false;
	if(s.equals("f")) return false;
	if(s.equals("n")) return false;
	if(s.equals("0")) return false;
	def msg =  "Doesn't look like a boolean : ["+ s + "]";
	log.warn(msg);
	throw new IllegalArgumentException(msg);
	}


def getKeyValue(meta,key,defValue) {
	if(!meta.containsKey(key)) {
		log.warn1("meta doesn't contains key \""+key+"\". Using default:\""+defValue+"\".",firstOnly: true);
		}
	return meta.getOrDefault(key,defValue);
	}

String getModules(String s) {
	String ret="";
	for(String s1:s.split("[ \t]+")) {
		if(s1.isEmpty()) continue;
		if(s1.contains("/")) {
			ret += " "+s1;
			continue;
			}
		if(s1.equals("plink")) {
			ret += " plink/1.90b6.18";
			continue;
			}
		if(s1.equals("gatk4")) {
			ret += " gatk/0.0.0";
			continue;
			}
		if(workflow.userName.equals("lindenbaum-p")) {
			if(s1.equals("picard") || s1.equals("bcftools") || s1.equals("samtools") || s1.equals("bedtools")) {
				ret += " "+s1+"/0.0.0";
				continue;
				}
			}
		ret+=" "+s1;
		}
	return ret;
	}

String moduleLoad(String s) {
	if(isBlank(s)) return "";
	return "module load " + getModules(s);
	}

boolean isHs37d5(String reference) {
	String fname = file(reference).getSimpleName().toLowerCase();
	if(fname.contains("hs37d5")) return true;
	return false;
	}

boolean isHg19(String reference) {
	if(isHs37d5(reference)) return true;
	String fname = file(reference).getSimpleName().toLowerCase();
	if(fname.contains("human_g1k_v37")) return true;
	if(fname.contains("hg19")) return true;
	if(fname.contains("grch37")) return true;
	return false;
	}

boolean isHg38(String reference) {
	String fname = file(reference).getSimpleName().toLowerCase();
	if(fname.contains("hs38me")) return true;
	if(fname.contains("hg38")) return true;
	if(fname.contains("grch38")) return true;
	return false;
	}


boolean isBlank(String s) {
	return s==null || s.trim().isEmpty() || s.equals(".");
	}

boolean isUrl(Object o) {
	if(o==null) return false;
	if(o instanceof java.net.URL) return true;
	final String s = String.valueOf(o);
	if(s.startsWith("https://") || s.startsWith("http://") || s.startsWith("ftp://")) return true;
	return false;
	}

String assertKnownReference(String reference) {
	if(isHg19(reference)) return reference;
	if(isHg38(reference)) return reference;
	throw new IllegalArgumentException("reference '${reference}' is unknown (or I cannot guess it from it's name). Currently supported references are defined in 'modules/utils/functions.nf'");
	}

String getUcscBuildName(String reference) {
	if(isHg19(reference)) return "hg19";
	if(isHg38(reference)) return "hg38";
	return "undefined";
	}

String getGencodeGtfUrl(String reference) {
	if(isHg19(reference)) return "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh37_mapping/gencode.v40lift37.annotation.gtf.gz";
	if(isHg38(reference)) return "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz";
	return "undefined";
	}

String getGencodeGff3Url(String reference) {
	if(isHg19(reference)) return "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh37_mapping/gencode.v40lift37.annotation.gff3.gz";
	if(isHg38(reference)) return "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gff3.gz";
	return "undefined";
	}

void assertTrue(boolean b,String msg) {
	if(!b) throw new IllegalStateException("assertion failed."+msg);
	}

void assertFalse(boolean b,String msg) {
	assertTrue(!b,msg);
	}

String assertNotEmpty(String s,String msg) {
	if(s==null || s.trim().isEmpty()) throw new IllegalStateException("Null or empty parameter:"+(msg==null?"a param is null":msg));
	return s;
	}

def assertFileExists(def f,String msg) {
	final String fstr = f==null?null:f.toString();
	if(isBlank(fstr)) throw new IllegalStateException("Missing file ("+fstr+"):"+msg==null?"a param is null":msg);
	final java.io.File file = new java.io.File(fstr);
	if(!file.exists()) throw new IllegalStateException("Missing File ("+fstr+"):"+(msg==null?"file is missing.":msg));
	return f;
	}


boolean hasFeature(Map meta,String v) {
        def disableFeatures = getKeyValue(meta,"disableFeatures","");
        if((disableFeatures instanceof java.lang.Boolean)) return true;
        final String[] tokens = disableFeatures.split("[, ;]+");
        for(String t: tokens) {
                t= t.trim();
                if(t.equalsIgnoreCase(v)) {
                        log.warn(v+" is disabled.");
                        return false;
                        }
                }
        return true;
        }

String getClassTaxonomy(def f) {
	if(f==null) return "null";
	String s="";
	def o = f.getClass();
	if(o==null) return "null";
	do {
		s+= o.getName()+";";
		o = o.getSuperclass();
		} while(o!=null && !(o.equals(java.lang.Object.class)));
	return s;
	}

String toAbsolutePath(def f) {
	if(f==null) throw new IllegalArgumentException("null path");
	if((f instanceof org.codehaus.groovy.runtime.GStringImpl)) {
		f = f.toString()
		}
	if((f instanceof String)) {
		if(!f.startsWith(java.io.File.separator)) throw new IllegalArgumentException("Not an absolute path:\""+f+"\"");
		return f;
		}
	if((f instanceof java.io.File)) return java.io.File.class.cast(f).getAbsolutePath();
	if((f instanceof java.nio.file.Path)) return java.nio.file.Path.class.cast(f).toAbsolutePath().toString();
	
	throw new IllegalArgumentException("Cannot convert "+f+" to path. ("+getClassTaxonomy(f)+")");
	}


String removeCommonSuffixes(String s) {
	List<String> L = [
		".gz",".txt",".csv",".idx",".tbi",".list",".tmp",".bgz",".zip",".tar",
		".bed",".interval_list",
		".vcf",".bcf",
		".fa",".fasta",".fna",
		".bam",".cram",".bai",".csi",".crai",
		".gtf",".gff",".gff3",
		".fq",".fastq",".ora"
		]
	boolean done = false;
	while(!done) {
		done = true;
		for(String x:L) {
			if(x.length()< s.length() && s.toLowerCase().endsWith(x)) {
				done=false;
				s=s.substring(0,s.length()-x.length());
				}
			}
		}
	return s;
	}


java.util.Map<String,Object> readContigFile(def f) {
	final java.util.Map<String,Object> hash = new HashMap<>();
	java.io.File fin = new java.io.File(toAbsolutePath(f));
	try(FileReader fr =new FileReader(fin)) {
		java.util.Properties props = new java.util.Properties();
		props.load(fr);
		hash.putAll(props);
		}
	catch(Throwable err) {
		log.warn(err.getMessage());
		err.printStackTrace();
		throw err;
		}
	return hash;
	}

String jvarkit(s) {
        return "\${JVARKIT_DIST}/"+s+".jar";
        }

String escapeXml(String s) {
	final StringBuilder sb= new StringBuilder(s.length());
	for(int i=0;i< s.length();i++) {
		char c = s.charAt(i);
		switch(c) {
			case '<': sb.append("&lt;");break;
			case '>': sb.append("&gt;");break;
			case '"': sb.append("&quot;");break;
			case '\'': sb.append("&apos;");break;
			case '&': sb.append("&amp;");break;
			default: sb.append(c);break;
			}
		}
	return sb.toString();
	}

String __getVersionCmd(java.util.Set<String> tools) {
	final StringBuilder sb= new StringBuilder("<dl>");
	for(String t0 : tools) {
		if(t0.isEmpty()) continue;
		sb.append("<dt>").append(t0).append("</dt><dd>");
		final String t = t0.toLowerCase();
		if(t.equals("bwa")) {
			sb.append("\$(bwa 2>&1 | grep Version)");
			}
		else if(t.equals("samtools")) {
			sb.append("\$(samtools  --version | head -n 1| cut -d ' ' -f2)");
			}
		else if(t.equals("bcftools")) {
			sb.append("\$(bcftools --version-only)");
			}
		else if(t.equals("bedtools")) {
			sb.append("\$(bedtools --version)");
			}
		else if(t.equals("gatk")) {
			sb.append("\$(gatk --version 2> /dev/null  | paste -s -d ' ' )");
			}
		else if(t.equals("snpeff")) {
			sb.append("\$(java -jar \${SNPEFF_JAR} -version)");
			}
		else if(t.equals("java")) {
			sb.append("\$( java -version 2>&1 | paste -s  -d ' ' )");
			}
		else if(t.equals("javac")) {
			sb.append("\$(javac -version 2>&1 )");
			}
		else if(t.equals("wget") || t.equals("git") || t.equals("tabix") || t.equals("bgzip") || t.equals("rvtest") || t.equals("gcc") || t.equals("g++") || t.equals("datamash")) {
			sb.append("\$("+t+" --version |head -n1)");
			}
		else if(t.equals("awk")) {
			sb.append("\$(awk --version | head -n1)");
			}
		else if(t.equals("cnvnator")) {
			sb.append("\$("+t+" 2>&1 | grep CNVn -m1)");
			}
		else if(t.equals("xsltproc")) {
			sb.append("\$("+t+" --version | paste -s -d ' ')");
			}
		else if(t.equals("gs")) {
			sb.append("\$(gs --version)");
			}
		else if(t.equals("r")) {
			sb.append("\$(R --version | head -n1)");
			}
		else if(t.equals("vep")) {
			sb.append("\$(vep --help 2>&1 | grep \"^  ensembl\" | tr -s \" \" | paste -s -d ' ')");
			}
		else if(t.equals("delly") || t.equals("delly2")) {
			sb.append("\$(delly --version  | head -n1 | cut -d':' -f 2)");
			}
		else if(t.startsWith("jvarkit/")) {
			final String j = t.substring(8);
			sb.append("\$(java -jar \${JVARKIT_DIST}/" + j + ".jar --version )");
			}
		else if(t.startsWith("picard/")) {
			final String j = t.substring(7);
			sb.append("\$(java -jar \${PICARD_JAR} "+t0+" --version 2>&1)");
			}
		else	{
			log.warn("getVersionCmd is undefined for  :"+t0);
			}
		sb.append("</dd>");
		}
	sb.append("</dl>");
	return sb.toString();
	}


String getVersionCmd(String s) {
	final java.util.Set<String> set = new java.util.TreeSet<>();
	boolean java=false;
	for(String t : s.trim().split("[ \t,;]+")) {
		set.add(t);
		if(t.startsWith("jvarkit/") || t.startsWith("picard/") || t.startsWith("gatk") || t.equalsIgnoreCase("snpeff")) {
			java=true;
			}
		}
	if(java) set.add("java");
	set.remove("");
	return __getVersionCmd(set);
	}

/** extract md5 cheksum from string */
def md5(def in) {
	final java.security.MessageDigest _md5;
	try {
		_md5 = java.security.MessageDigest.getInstance("MD5");
	} catch (final java.security.NoSuchAlgorithmException e) {
		throw new RuntimeException("MD5 algorithm not found", e);
		}
	
	_md5.reset();
	_md5.update(in.toString().getBytes());
	String s = new java.math.BigInteger(1, _md5.digest()).toString(16);
	if (s.length() != 32) {
		final String zeros = "00000000000000000000000000000000";
		s = zeros.substring(0, 32 - s.length()) + s;
		}
	return s;
	}


void _dumpParams(final int level,final String margin, final java.util.Map map, java.lang.Appendable out) {
	 def fmt = org.fusesource.jansi.Ansi.ansi()
	for(Object k: map.keySet()) {
		final Object v = map.get(k);
		out.append(margin);
		out.append(" + ");

		//fmt = fmt.fg(Color.GREEN);
		out.append(k.toString());

		//fmt = fmt.fg(Color.DEFAULT)	

		out.append(" : ");
		if(v instanceof Map) {
			if(level==0 && k.equals("genomes")) {
				out.append(" (hidden) ....\n");
				}
			else {
				out.append("\n");
				_dumpParams(level+1, margin+"    ",(Map)v, out);
				}
			}
		else
			{
			out.append(String.valueOf(v));
			out.append("\n");
			}
		}
	}

void dumpParams(final java.util.Map map) {
	System.out.println("# Params:");
	_dumpParams(0,"",map,System.out);
	System.out.println();
	}

String paramsToString(final java.util.Map map) {
	StringBuilder sb = new StringBuilder();
	_dumpParams(0,"",map,sb);
	sb.append("\n");
	return sb.toString();
	}


Object slurpJsonFile(f) {
	if(f instanceof nextflow.processor.TaskPath) {
		return slurpJsonFile(f.toRealPath())
		}
	def slurper = new groovy.json.JsonSlurper();
	return slurper.parse(f)
	}


Map assertKeyExists(final Map hash,final String key) {
    if(hash[key]==null) throw new IllegalArgumentException("no key '${key}' in ${hash}.");
    return hash;
}

boolean testKeyExistsAndNotEmpty(final Map hash,final String key) {
  def value = hash.get[key];
  if(value==null || value.isEmpty() || value.equals(".")) return false;
  return true;
}

Map assertKeyExistsAndNotEmpty(final Map hash,final String key) {
    assertKeyExists(hash,key);
    def value = hash[key];
    if(value==null || value.isEmpty()) throw new IllegalArgumentException("empty ${key}'in ${hash}");
    return hash;
}

Map assertKeyMatchRegex(final Map hash,final String key,final String regex) {
    assertKeyExists(hash,key);
    def value = hash.get(key);
    if(!value.matches(regex)) throw new IllegalArgumentException(" ${key}'in ${hash} doesn't match regex '${regex}'.");
    return hash;
  }

/** return true if error looks like at error from the cluster, not from the script itself */
boolean isClusterError(task) {
	if(task==null) return false;
	if(task.previousException == null) return false;
	String s = task.previousException.message;
	if(s.startsWith("Error submitting process ") && s.endsWith("for execution")) {
		return true;
		}
	return false;
	}

/** called by 'makeKey' below, private use  */
void _makeKey(digest,o) {   
    if(o==null) {
        digest.update((byte)0);
        }
    else if(o instanceof List) {
        digest.update((byte)'[');
        for(int i=0;i< o.size();i++) {
            _makeKey(digest,o[i]);
            if(i>0) digest.update((byte)',');
        }
       digest.update((byte)']');
    } else if(o instanceof Map) {
        def keys = o.keySet().sort()
        digest.update((byte)'[');
        for(int i=0;i< keys.size();i++) {
            _makeKey(digest,keys[i]);
             digest.update((byte)':');
             _makeKey(digest,o.get(keys[i]));
            if(i>0) digest.update((byte)',');
        }
       digest.update((byte)']');
    } else if(o instanceof java.lang.String) {
            digest.update(o.getBytes("UTF-8"));
    } else if(o instanceof java.lang.Number || o instanceof java.lang.Boolean) {
		_makeKey(digest,o.getClass().toString());
        _makeKey(digest,o.toString());
    } else if(o instanceof java.nio.file.Path) {
		   _makeKey(digest,o.getClass().toString());
           if(workflow.stubRun && !o.exists())
		{
		 _makeKey(digest,o.toString());
		}
	    else
		{
	           _makeKey(digest,o.toRealPath().toString());
		}
    } else {
       throw new IllegalArgumentException("cannot invoke _makeKey with class=${o.getClass()}");
    }
}

/** convert a java/groovy object to a md5 checksum */
String makeKey(o) {
    java.security.MessageDigest md5 = java.security.MessageDigest.getInstance("MD5");
    _makeKey(md5,o);
    java.math.BigInteger bigInt = new java.math.BigInteger(1, md5.digest());
    String s= "K"+bigInt.toString(16);
    return s;
}


/** test if a file contains a data after uncompressing if needed */
boolean isEmptyGz(def p) {
        if(workflow!=null && workflow.stubRun==true) return false;
	verify(p!=null,"isEmptyGZ: p is null?");
        if(!(p instanceof java.nio.file.Path)) {
            throw new IllegalArgumentException("expected a path for ${p} but got ${p.class}");
            }
        if(!java.nio.file.Files.exists(p)) {
            throw new java.io.FileNotFoundException("not found ${p}");
            }
        if(p.getFileName().toString().toLowerCase().endsWith(".gz")) {
            try (java.util.zip.GZIPInputStream in = new  java.util.zip.GZIPInputStream(java.nio.file.Files.newInputStream(p)) ) {
                return in.read()==-1;
                }
            }
        else
            {
            try (java.io.InputStream in = java.nio.file.Files.newInputStream(p)) {
                return in.read()==-1;
                }
            }
        }


/** flat map by array index */
List flatMapByIndex(L0,index) {
	verify(L0!=null,"flatMapByIndex L0 is null");
	verify((L0 instanceof List),"flatMapByIndex should be a List ${L0.class}");
	verify(index>=0 && index<L0.size() , "flatMapByIndex index=0<=${index}<${L0.size()}. out or range");
	def pivot = L0[index];
	verify((pivot instanceof List), "flatMapByIndex index in ${L0}. pivot is NOT a list ${pivot.class}");

	try {
	def L1 = [];
	for(int i=0;i< pivot.size();i++) {
		def L2=[];
		for(int j=0;j< L0.size();j++) {
			def value = null;
			if(j==index) {
				value = pivot[i];
				}
			else
				{
				value = L0[j];
				}
			L2.add(value==null?null:(value instanceof java.lang.Cloneable?value.clone():value));
			}
		L1.add(L2);
		}
	return L1;
	} catch(Throwable err) {
		log.warn("flatMapByIndex: ${err}");
		throw err;
		}
	}

/** split file/index using the htslib method . If found return array with two elements (file/index) or one element (file) */
List htslibSplitIndex(String fname) {
	String delim = "##idx##";
	int i = fname.indexOf(delim);
	if(i<1) return [fname];
	String f1 = fname.substring(0,i);
	String f2 = fname.substring(i+delim.length());
	return [f1,f2];
}


/** extract sample from illumina fastq name : e.g. S1_S10_L008_R2_001.fastq.ora -> S1 */
Map extractIlluminaName(String f) {
	try {
		f = file(f).name;
        // _S\\d_L\\d+_R\\d_\\d+
		java.util.regex.Pattern pattern = java.util.regex.Pattern.compile("^(.*)_(S\\d+)_(L\\d+)_(R\\d)_(\\d+)\\.(fq|fastq)(\\.(ora|gz))?");
		java.util.regex.Matcher matcher = pattern.matcher(f);
		if(!matcher.find()) return null;
		
		String sn = matcher.group(1);
		verify(sn.indexOf("/")==-1,"Whattt ? ${f}");
		return [
			id: sn,
			sample: sn,
			index : matcher.group(2),
			lane: matcher.group(3),
			side : matcher.group(4),
			split :  matcher.group(5)
			];
            
		}
	catch(Throwable err) {
		err.printStackTrace();
		log.warn("cannot extractORASampleName "+f+" "+err.getMessage());
		return null;
		}
	}

/** assertion checker */
boolean verify(boolean t, String msg) {
	if(t) return true;
	log.warn(String.valueOf(msg));
	throw new IllegalStateException(String.valueOf(msg));
	}

/** media https://stackoverflow.com/a/44376148/58082 */
def median(data) {
    def copy = data.toSorted()
    def middle = data.size().intdiv(2)

    // you can omit the return in groovy for the last statement
    data.size() % 2 ? copy[middle] : (copy[middle - 1] + copy[middle]) / 2
   }


void runOnComplete(def wf) {
wf.onComplete {
    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${wf.complete}
        Duration    : ${wf.duration}
        Success     : ${wf.success}
        workDir     : ${wf.workDir}
        exit status : ${wf.exitStatus}
        """ : """
        Failed: ${wf.errorReport}
        exit status : ${wf.exitStatus}
        """
    )
}
}
