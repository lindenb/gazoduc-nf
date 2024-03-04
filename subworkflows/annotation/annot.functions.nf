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
        def with = "with_"+key
        if(!params.containsKey("annotations")) return false;
        if(!params.annotations.containsKey(with)) {
                log.warn("params.annotations."+with+" missing. return false;");
                return false;
                }
        if(isBlank(params.annotations,with)) return false;
        return (params.annotations[with] as boolean)
        }

boolean isSoftFilter(key) {
        def f = "softfilter_"+key
        if(!params.containsKey("annotations")) return true;
        if(!params.annotations.containsKey(f)) {
                log.warn("params.annotations."+f+" missing. return true;");
                return true;
                }
        if(isBlank(params.annotations,f)) return true;
        return (params.annotations[f] as boolean)
        }

String backDelete(json) {
	return "## rm -f \"${json.vcf}\"  \"${json.index}\" ";
	}