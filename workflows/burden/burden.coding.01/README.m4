m4_define(`MACRO_WORKFLOW',burdencoding)m4_dnl
m4_define(`MACRO_MAIN_NF',gazoduc-nf/burden/burden.coding.01/burden.coding.01.nf)m4_dnl
m4_define(`MACRO_MAIN_CFG',../../../confs/by_workflow/burden.coding.01.cfg)m4_dnl
m4_sinclude(../../../src/main/m4/common.m4)m4_dnl
MACRO_H1(MACRO_WORKFLOW)

MACRO_H2(About)

Burden test for rare variants in coding regions.

MACRO_HOW_TO_INSTALL

MACRO_PEDIGREE

MACRO_H2(Parameters)

MACRO_PARSE_CONFIG(`../../../lib/nfconfigparser.jar',` ../../../confs/by_workflow/wgselect.basic.cfg')

MACRO_EXECUTE(	--vcf /path/to/vcf \
	MACRO_COMMON_PARAMS \
	--pedigree "/path/to/file.ped" \
	--bed /path/to/input.bed 
)

MACRO_FOOTER
