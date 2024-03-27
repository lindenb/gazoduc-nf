
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
import java.io.*;
import java.util.*;
import java.util.stream.*;

import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import com.github.lindenb.jvarkit.iterator.EqualIterator;
import com.github.lindenb.jvarkit.lang.StringUtils;

public class Minikit {
private static final boolean WITH_BND = __WITH_BND__;
private final Map<String,String> convertHash = new HashMap<>();
private static final int MIN_SV_LEN = __MIN_SV_LEN__;
private static final int MAX_SV_LEN = __MAX_SV_LEN__;

 
private String convert(final String s) {
	return convertHash.get(s);
	}
private String type2color(final String s) {
	if(s.equals("DEL")) return "0,0,255";
	if(s.equals("INS")) return "0,255,0";
	if(s.equals("INV")) return "255,0,0";
	if(s.equals("BND")) return "255,255,0";
	if(s.contains("DUP")) return "255,0,255";
	return "0,0,0";
	}
// A]hs37d5:12060965]]
private List<Locatable> parseBnd(final VariantContext ctx) {
	if(!WITH_BND) return Collections.emptyList();
	return ctx.getAlternateAlleles().stream().
		filter(A->A.isSymbolic()).
		map(A->A.getDisplayString()).
		flatMap(S->Arrays.stream(S.split("[,]"))).
		map(S->{
		int x1 = S.indexOf("[");
		if(x1==-1) x1= S.indexOf("]");
		if(x1==-1) return null;
		int x2 =  S.indexOf("[",x1+1);
		if(x2==-1) x2 =S.indexOf("]",x1+1);
		if(x2==-1) return null;
		int colon = S.indexOf(":",x1+1);
		if(colon==-1 || colon>=x2) return null;
		final String contig = convert(S.substring(x1+1,colon));
		if(contig==null) return null;
		int pos = Integer.parseInt(S.substring(colon+1,x2));
		return new Interval(contig,pos,pos);
		}).
		filter(R->R!=null).
		collect(Collectors.toList());
	}

private int getSVLen(final VariantContext ctx) {
    final int svlen;
	if(ctx.hasAttribute("SVLEN")) {
		svlen = Math.abs( ctx.getAttributeAsInt("SVLEN",0));
		}
	else if(ctx.hasAttribute("SVINSLEN")) {
		svlen = Math.abs( ctx.getAttributeAsInt("SVINSLEN",0));
		}
	else if(ctx.hasAttribute("SVINSSEQ")) {
		svlen = Math.abs( ctx.getAttributeAsString("SVINSSEQ","").length());
		}
	else
		{
		svlen = ctx.getLengthOnReference();
		}
    return svlen;
    }

private String getSVType(VariantContext ctx) {
	return  ctx.getAttributeAsString("SVTYPE",".");
	}

private List<List<VariantContext>>  groupSame(List<VariantContext> L)  {
    if(L.size()<2) return Collections.singletonList(L);
    final  List<List<VariantContext>> ret = new ArrayList<>();
    L = new ArrayList<>(L);
    while(!L.isEmpty()) {
        List<VariantContext> L2 = new ArrayList<>();
        ret.add(L2);
        final VariantContext curr = L.remove(0);
        L2.add(curr);
        int i=0;
        while(i< L.size()) {
            final VariantContext item = L.get(i);
            if(item.contigsMatch(curr) &&
               curr.getStart()==item.getStart() &&
               getSVLen(curr)==getSVLen(item) &&
               getSVType( curr).equals(getSVType(item)) &&
               (!getSVType(curr).equals("BND") || parseBnd(curr).equals(parseBnd(item)))
               ) {
                L2.add(L.remove(i));
                }
            else {
                i++;
                }
            }
        }
    return ret;
    }

private int doWork(final List<String> args) {
	if(args.size()!=3) {
		System.err.println("usage vcfname out.bed out.bedpe");
		System.exit(-1);
		}
	try
		{
		for(int i=1;i<=22;i++) {
			convertHash.put(""+i,"chr"+i);
			convertHash.put("chr"+i,"chr"+i);
			}

		convertHash.put("X","chrX");
		convertHash.put("chrX","chrX");
		convertHash.put("Y","chrY");
		convertHash.put("chrY","chrY");
		convertHash.put("MT","chrM");
		convertHash.put("M","chrM");

		String vcfName= new File(args.get(0)).getName();
		if(vcfName.endsWith(".vcf.gz")) vcfName=vcfName.substring(0, vcfName.length() -7);
		if(vcfName.endsWith(".vcf")) vcfName=vcfName.substring(0, vcfName.length() -4);
		if(vcfName.endsWith(".bcf")) vcfName=vcfName.substring(0, vcfName.length() -4);

		try(PrintWriter pw1 = new PrintWriter(args.get(1)); PrintWriter pw2 = new PrintWriter(args.get(2))){
			try(final VCFIterator r0 = new VCFIteratorBuilder().open(System.in)) {
				final List<String> samples = r0.getHeader().getSampleNamesInOrder();
				final SAMSequenceDictionary dict = r0.getHeader().getSequenceDictionary();
				final EqualIterator<VariantContext> r = new EqualIterator<>(r0,(A,B)->{
				       int i= A.getContig().compareTo(B.getContig());
				       if(i!=0) return i;
				       return Integer.compare(A.getStart(),B.getStart()); 
				       });
				
				while(r.hasNext()) {
				for(List<VariantContext> array : groupSame(r.next())) {
					final VariantContext ctx = array.get(0);
					final SAMSequenceRecord ssr = (dict==null?null:dict.getSequence(ctx.getContig()));
					if(ssr==null) continue;
					final int cipos = array.stream().flatMap(V->V.getAttributeAsIntList("CIPOS",0).stream()).mapToInt(V->V.intValue()).min().orElse(0);
					final int ciend = array.stream().flatMap(V->V.getAttributeAsIntList("CIEND",0).stream()).mapToInt(V->V.intValue()).max().orElse(0);

					int start0 = Math.max(0,(ctx.getStart()-1) + cipos);
					int end0 = ctx.getEnd() + ciend;

					if(start0>=end0 || end0>ssr.getSequenceLength()) {
						start0 = ctx.getStart()-1;
						end0 = ctx.getEnd();
						}

					final String type= getSVType(ctx);
					if(StringUtils.isBlank(type)) continue;
					String contig =  convert(ctx.getContig());
					if(StringUtils.isBlank(contig)) continue;
					final Set<String> samplesSet = array.stream().flatMap(G->G.getGenotypes().stream()).
							filter(G->G.hasAltAllele()).
							map(G->G.getSampleName()).
							collect(Collectors.toCollection(TreeSet::new));
					
					
					if(type.equals("BND")) {
						for(Locatable bnd : parseBnd(ctx)) {
							final String contig2 =  convert(bnd.getContig());
							if(contig2==null || contig2.isEmpty()) continue;
							pw2.print(contig);
							pw2.print("\t");
							pw2.print(start0);
							pw2.print("\t");
							pw2.print(end0);
							pw2.print("\t");
							pw2.print(samplesSet.isEmpty()?".":(samplesSet.size()>20?"N="+samplesSet.size():String.join("|", samplesSet)));
							pw2.print("\t");
							pw2.print(1000-(int)(1000.0*(samplesSet.size()/(double)samples.size())));
							pw2.print("\t");
							pw2.print(array.stream().allMatch(V->V.isFiltered())?"1":"2");//value
							pw2.print("\t");
							pw2.print(".");//exp
							pw2.print("\t");
							pw2.print(contig.equals("contig2")?"0,0,255":"255,0,0");//color
							pw2.print("\t");
							pw2.print(contig);//source chrom
							pw2.print("\t");
							pw2.print(ctx.getStart()-1);//source start
							pw2.print("\t");
							pw2.print(ctx.getEnd());//source end
							pw2.print("\t");
							pw2.print(".");//source name
							pw2.print("\t");
							pw2.print(".");//source strand
							pw2.print("\t");
							pw2.print(contig2);//target chrom
							pw2.print("\t");
							pw2.print(bnd.getStart()-1); //target start
							pw2.print("\t");
							pw2.print(bnd.getEnd());//target end
							pw2.print("\t");
							pw2.print(".");//target name
							pw2.print("\t");
							pw2.print(".");//target strand
							pw2.print("\t");
							pw2.print(StringUtils.ifBlank(array.stream().filter(V->V.isFiltered()).flatMap(V->V.getFilters().stream()).collect(Collectors.toSet()).stream().collect(Collectors.joining(",")),"PASS"));//filters
							pw2.print("\t");
							pw2.print(StringUtils.ifBlank(array.stream().flatMap(G->G.getGenotypes().stream()).filter(G->G.hasAltAllele()).map(t->t.getType().name()).collect(Collectors.toSet()).stream().collect(Collectors.joining(",")),"."));//gt type
							pw2.print("\t");
							pw2.print(vcfName);//vcf Name
							pw2.println();
							}//end of BND
						} else //not BND
						{
						final int svlen = getSVLen(ctx);
						
						if( svlen < MIN_SV_LEN) continue;
						if( svlen > MAX_SV_LEN) continue;

						pw1.print(contig);
						pw1.print("\t");
						pw1.print(start0);
						pw1.print("\t");
						pw1.print(end0);
						pw1.print("\t");
						pw1.print(samplesSet.isEmpty()?".":(samplesSet.size()>20?"N="+samplesSet.size():String.join("|", samplesSet)));
						pw1.print("\t");
						pw1.print(1000-(int)(1000.0*(samplesSet.size()/(double)samples.size())));
						pw1.print("\t");
						pw1.print("+");
						pw1.print("\t");
						pw1.print(ctx.getStart()-1);
						pw1.print("\t");
						pw1.print(ctx.getEnd());
						pw1.print("\t");
						pw1.print(type2color(type));
						pw1.print("\t");
						pw1.print(type);//type
						pw1.print("\t");
						pw1.print(svlen);//svlen
						pw1.print("\t");
						pw1.print(StringUtils.ifBlank(array.stream().filter(V->V.isFiltered()).flatMap(V->V.getFilters().stream()).collect(Collectors.toSet()).stream().collect(Collectors.joining(",")),"PASS"));//filters
						pw1.print("\t");
						pw1.print(StringUtils.ifBlank(array.stream().flatMap(G->G.getGenotypes().stream()).filter(G->G.hasAltAllele()).map(t->t.getType().name()).collect(Collectors.toSet()).stream().collect(Collectors.joining(",")),"."));//gt type
						pw1.print("\t");
						pw1.print(vcfName);//vcf Name
						pw1.println();
						}
					} // end groupSame
				}//end r.hasNext
		    r.close();
			} // end r0
		pw1.flush();
		pw2.flush();
		}
	catch(final Throwable err) {
		err.printStackTrace();
		return -1;
		}

	return 0;
	}
catch(final Throwable err2) {
	err2.printStackTrace();
	return -1;
	}
}

public static void main(final String args[]) {
	try {
		int ret= new Minikit().doWork(Arrays.asList(args));
		System.exit(ret);
		}
	catch(Throwable err) {
		err.printStackTrace();
		System.exit(-1);
		}
	}

}
