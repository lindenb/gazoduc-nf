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
	def o = meta.get(key);
	if(o == null) return false;
	if((o instanceof Boolean)) {
		return (Boolean)o;
		}
	final String s= String.valueOf(o).trim();
	if(s.equalsIgnoreCase("true")) return true;
	if(s.equalsIgnoreCase("yes")) return true;
	if(s.equalsIgnoreCase("T")) return true;
	if(s.equalsIgnoreCase("y")) return true;
	if(s.equalsIgnoreCase("1")) return true;

	if(s.equalsIgnoreCase("false")) return false;
	if(s.equalsIgnoreCase("no")) return false;
	if(s.equalsIgnoreCase("F")) return false;
	if(s.equalsIgnoreCase("n")) return false;
	if(s.equalsIgnoreCase("0")) return false;

	//if(__warnings.add(key)) {
		log.warn("meta doesn't contains boolean key \""+key+"\". Using default:\"false\".");
		//}
	return false;
	}


def getKeyValue(meta,key,defValue) {
	if(!meta.containsKey(key)) {
		log.warn("meta doesn't contains key \""+key+"\". Using default:\""+defValue+"\".");
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

boolean isHg19(String reference) {
	String fname = file(reference).getSimpleName();
	if(fname.contains("hs37d5")) return true;
	if(fname.contains("human_g1k_v37")) return true;
	if(fname.contains("hg19")) return true;
	return false;
	}

boolean isHg38(String reference) {
	String fname = file(reference).getSimpleName();
	if(fname.contains("hs38me")) return true;
	return false;
	}

