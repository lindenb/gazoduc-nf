
private boolean acceptControl(final VariantContext vc,String sm) {
	final Genotype g = vc.getGenotype(sm);
	if(g!=null && g.isHomVar()) return false;
	return true;
	}



private boolean acceptTrio(final VariantContext vc,String cm,String fm,String mm) {
	final Genotype c = vc.getGenotype(cm);
    if(c==null) return false;
	final Genotype m = vc.getGenotype(mm);
	final Genotype f = vc.getGenotype(fm);
	// child must be HET
	if(!c.isHet()) return false;
	// discard de novo
	if(f.isHomRef() && m.isHomRef()) return false;
	// parent shouldn't be homvar
	if(f.isHomVar()) return false;
	if(m.isHomVar()) return false;
	// at least one parent should be het
	if(!(f.isHet() || m.isHet())) return false;
	return true;
	}

public Object apply(final VariantContext variant) {
    m4_include(custom.m4)
    return false;
    }
