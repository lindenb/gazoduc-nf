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
	try(BufferedReader br=java.nio.file.Files.newBufferedReader(fai)) {
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



boolean test_fai_chr1(fai,K1_length) {
	try(BufferedReader br=  java.nio.file.Files.newBufferedReader(fai)) {
		for(;;) {
			final String line = br.readLine();
			if(line==null) break;
			final String[] tokens =line.split("[\t]");
			if(tokens[0].equals("chr1") || tokens[0].equals("1")) {
				return tokens[1].equals(K1_length);
			}
		}
	} catch(Throwable err) {
		log.info("Error err: "+err);
		err.printStackTrace();
		throw err;
	}
	return false;
}

boolean isGRCh38(fai) {
	def k1 = k1_signature();
	return test_fai_chr1(fai,String.valueOf(k1.hg38));
}

boolean isGRCh37(fai) {
	def k1 = k1_signature();
	return test_fai_chr1(fai,String.valueOf(k1.hg19));
}