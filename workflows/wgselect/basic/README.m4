m4_define(`MACRO_WORKFLOW',wgselect)m4_dnl
m4_define(`MACRO_MAIN_NF',gazoduc-nf/workflows/wgselect/basic/wgselect.basic.nf)m4_dnl
m4_define(`MACRO_MAIN_CFG',../../../confs/by_workflow/wgselect.basic.cfg)m4_dnl
m4_sinclude(../../../src/main/m4/common.m4)m4_dnl
MACRO_H1(MACRO_WORKFLOW)

MACRO_H2(About)

This workflow filter out bad quality variants from a VCF file.

MACRO_HOW_TO_INSTALL

MACRO_PEDIGREE

MACRO_H2(Parameters)

MACRO_OPTIONS(`MACRO_AWK_PARAMS(../../../confs/by_workflow/wgselect.basic.cfg,)m4_dnl
MACRO_AWK_PARAMS(../../../confs/default.params.cfg,)m4_dnl
MACRO_AWK_PARAMS(../../../confs/genomeId.params.cfg,)m4_dnl
MACRO_AWK_PARAMS(../../../confs/by_subworkflow/wgselect.config,wgselect.)m4_dnl
')

MACRO_EXECUTE(	--vcf /path/to/vcf \
	MACRO_COMMON_PARAMS \
	--pedigree "/path/to/file.ped" \
	--bed /path/to/input.bed 
)

MACRO_FOOTER
