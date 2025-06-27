
process HET_COMPOSITE {
tag "${meta.id}"
label "process_single"
afterScript "rm -rf TMP"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta4),path(pedigree)
    tuple val(meta),path(vcf),path(vcfidx) // MUST BE ANNOTATED WITH SNPEFF, REMOVE FREQUENT VARIANTS
output:
	tuple val(meta),path("*.bcf"),path("*.bcf.csi"),emit:vcf
	tuple val(meta),path("genes.report"),emit:genes_report
	tuple val(meta),path("variants.report"),emit:variants_report
script:
    def prefix = task.ext.prefix?:vcf.simpleName+".hetcomposite"
"""
hostname 1>&2

mkdir -p TMP
set -x


cat << EOF > TMP/jeter.code

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
EOF

## all other samples are controls
comm -13 \
	<(cut -f 2 '${pedigree}' | sort | uniq) \
	<(bcftools query -l '${vcf}'| sort | uniq) |\
	awk '{printf("if(!acceptControl(variant,\\"%s\\")) return false;\\n",\$1);}' >> TMP/jeter.code

awk -F '\t' '(\$6=="control" || \$6=="unaffected") {printf("if(!acceptControl(variant,\\"%s\\")) return false;\\n",\$2);}'  '${pedigree}' >> TMP/jeter.code

awk -F '\t' '((\$6=="case" || \$6=="affected") && \$3!="0" && \$4!="0") {printf("if(acceptTrio(variant,\\"%s\\",\\"%s\\",\\"%s\\")) return true;\\n",\$2,\$3,\$4);}'  '${pedigree}'  >> TMP/jeter.code


echo "return false;}" >> TMP/jeter.code


bcftools view "${vcf}"  |\\
	jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP vcffilterjdk --body -f TMP/jeter.code > TMP/jeter2.vcf
mv TMP/jeter2.vcf TMP/jeter1.vcf




jvarkit -Xmx${task.memory.giga}G  -Djava.io.tmpdir=TMP vcfcomposite \\
	--extractors ANN/GeneId \\
	--filter "" \\
	--genes genes.report \\
	--pedigree  "${pedigree}" \\
	--report variants.report \\
	--tmpDir TMP \\
	--max-variants 30 \\
	TMP/jeter1.vcf > TMP/jeter2.vcf

mv TMP/jeter2.vcf TMP/jeter1.vcf


bcftools sort -T TMP/sort -o TMP/${prefix}.bcf -O b TMP/jeter1.vcf
bcftools index --threads ${task.cpus} TMP/${prefix}.bcf

mv TMP/${prefix}.bcf ./
mv TMP/${prefix}.bcf.csi ./

"""
}

