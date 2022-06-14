java.util.Set<String> __warnings = new java.util.HashSet<>();


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

	if(__warnings.add(key)) {
		log.warn("meta doesn't contains boolean key \""+key+"\". Using default:\"false\".");
		}
	return false;
	}


def getKeyValue(meta,key,defValue) {
	if(!meta.containsKey(key) && __warnings.add(key)) {
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

