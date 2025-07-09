boolean debug = false;
private boolean acceptControl(final VariantContext vc,String sm) {
	final Genotype g = vc.getGenotype(sm);
	if(g!=null && g.hasAltAllele()) {
        if(debug) System.err.println("control has alt "+sm);
        return false;
        }
	return true;
	}



private boolean acceptTrio(final VariantContext vc,String cm,String fm,String mm) {
	final Genotype c = vc.getGenotype(cm);
    if(c==null) return false;
	final Genotype m = vc.getGenotype(mm);
	final Genotype f = vc.getGenotype(fm);
	// child must have alt
	if(!c.hasAltAllele()) return false;
	// discard parent
    if(m!=null && m.hasAltAllele()) {
        if(debug) System.err.println("reject because father"+fm);
        return false;
        }
	if(f!=null && f.hasAltAllele()) {
        if(debug) System.err.println("reject because mother"+mm);
        return false;
        }
    if(debug) System.err.println("accept "+cm);
	return true;
	}

public Object apply(final VariantContext variant) {
    if(variant.hasAttribute("GNOMAD_AF")) {
        double af = variant.getAttributeAsDouble("GNOMAD_AF","")
        if(af >= 0.01) return false;
    }


    String svType = variant.getAttributeAsString("SVTYPE","");
    if(svType!=null && svType.equals("BND")) {
       if(debug) System.err.println("svType "+svType);
    }


    int svlen ;
    if(variant.hasAttribute("SVLEN")) {
        svlen = variant.getAttributeAsInt("SVLEN",0);
        }
    else
        {
        svlen = (variant.getEnd() - variant.getStart()) +1;
        }
    if(svlen<0) svlen = svlen*-1;

    if(svlen < __MIN_LEN__ ) {
        if(debug) System.err.println("small "+svlen);
        return false;
        }
    if(svlen > __MAX_LEN__ )  {
        if(debug) System.err.println("large "+svlen);
        return false;
        }

    final Set<String> children =new HashSet<>();

    m4_include(custom.m4)
    
    if(children.isEmpty()) {
        if(debug) System.err.println("no child found");
        return false;
        }
    try(java.io.PrintWriter w=new java.io.PrintWriter(new java.io.FileWriter("TMP/jeter.bed",true))) {
        for(String child: children) {
            w.println(variant.getContig()+"\t"+(variant.getStart()-1)+"\t"+variant.getEnd()+"\t"+svType+"\t"+svlen+"\t"+child);
            }
        w.flush();
        }
    catch(java.io.IOException err) {
        throw new RuntimeException(err);
        }

    return true;
    }