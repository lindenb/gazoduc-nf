/*

Copyright (c) 2024 Pierre Lindenbaum

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
boolean isBlank(hash,key) {
        if(hash==null) throw new IllegalArgumentException("[isBlank] hash is null");
        if(!hash.containsKey(key)) return true;
        def v  = hash.get(key);
        if(v==null || v.toString().trim().isEmpty()) return true;
        return false;
        }

boolean hasFeature(key) {
        if(!params.containsKey("annotations")) return false;
        if(!params.annotations.containsKey(key)) {
                log.warn("params.annotations."+key+" missing. return false;");
		return false;
		}
        if(!params.annotations[key].containsKey("enabled")) {
                log.warn("params.annotations."+key+".enabled missing. return false;");
                return false;
                }
        return (params.annotations[key].enabled.toBoolean());
        }

boolean isSoftFilter(key) {
	return !isHardFilter(key);
	}


boolean isHardFilter(key) {
        if(!params.containsKey("annotations")) return false;
        if(!params.annotations.containsKey(key)) {
                log.warn("params.annotations."+key+" missing. return false;");
		return false;
		}
        if(!params.annotations[key].containsKey("hard_filter")) {
                log.warn("params.annotations."+key+".hard_filter missing. return false;");
                return false;
                }
	
        return (params.annotations[key].hard_filter.toBoolean())
        }

String backDelete(json) {
	//return "## rm -f \"${json.vcf}\"  \"${json.index}\" ";
	return "\n### TODO\n";
	}

boolean isHg19(genomeId) {
	if(isBlank(params.genomes[genomeId],"ucsc_name")) return false;
	def u = params.genomes[genomeId].ucsc_name;
	return u.equals("hg19");
	}

boolean isHg38(genomeId) {
	if(isBlank(params.genomes[genomeId],"ucsc_name")) return false;
	def u = params.genomes[genomeId].ucsc_name;
	return u.equals("hg38");
	}

String hgName(genomeId) {
	if(isHg19(genomeId)) return "hg19";
	if(isHg38(genomeId)) return "hg38";
	return "";
	}
