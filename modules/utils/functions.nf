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


//static java.util.Set<String> __warnings = new java.util.HashSet<>();


/** return boolean from hasmap */
def getBoolean(meta,key) {
	if(!meta.containsKey(key)) return false;
	return parseBoolean(meta.get(key));
	}


def parseBoolean(def o) {
	if(o == null) return false;
	if((o instanceof Boolean)) {
		return (Boolean)o;
		}
	final String s= String.valueOf(o).trim();
	if(s.equalsIgnoreCase("true")) return true;
	if(s.equalsIgnoreCase("yes")) return true;
	if(s.equalsIgnoreCase("on")) return true;
	if(s.equalsIgnoreCase("T")) return true;
	if(s.equalsIgnoreCase("y")) return true;
	if(s.equalsIgnoreCase("1")) return true;

	if(s.equalsIgnoreCase("false")) return false;
	if(s.equalsIgnoreCase("no")) return false;
	if(s.equalsIgnoreCase("off")) return false;
	if(s.equalsIgnoreCase("F")) return false;
	if(s.equalsIgnoreCase("n")) return false;
	if(s.equalsIgnoreCase("0")) return false;

	log.warn("Doesn't look like a boolean : ["+ s + "]");
	return false;
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
			ret += "plink/1.90b6.18";
			continue;
			}
		if(s1.equals("gatk4")) {
			ret += "gatk/0.0.0";
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

boolean isHg19(String reference) {
	String fname = file(reference).getSimpleName().toLowerCase();
	if(fname.contains("hs37d5")) return true;
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

boolean isBlank(def s) {
	return s==null || s.toString().trim().isEmpty();
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

def getGnomadExomePath(def params,def reference) {
	if(params.containsKey("gnomad_exome_path")) {
		return params.get("gnomad_exome_path");
		}
	if(isHg19(reference)) return params.gnomad_exome_hg19;
	return "";
	}

def getGnomadGenomePath(def params,def reference) {
	if(params.containsKey("gnomad_genome_path")) {
		return params.get("gnomad_genome_path");
		}
	if(isHg19(reference)) return params.gnomad_genome_hg19;
	if(isHg38(reference)) return params.gnomad_genome_hg38;
	return "";
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

