/** signature of the chr1 of several build **/
def k1_signature() { return [
	hg19 : 249250621 ,
	hg38 : 248956422,
	canFam3 : 122678785,
	canFam4 : 123556469
	];
	}

def buildFromFai(fai) {
	def k1 = k1_signature();
	def hash = [:];
	for(key in k1.keySet()) {
		hash.put(k1[key].toString(),key);
		}
	try(BufferedReader br=Files.newBufferedReader(fai)) {
		String line;
		while((line=br.readLine())!=null) {
			String[] tokens = line.split("[\t]");
			if(!tokens[0].matches("(chr)?1")) continue;
			String build  = hash.get(tokens[1]);
			if(build!=null) return build;
		}
	}
	return "undefined_build";
}