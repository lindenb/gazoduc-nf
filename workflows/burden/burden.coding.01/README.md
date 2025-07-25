MACRO_H1(burdencoding)

MACRO_H2(About)

Burden test for rare variants in coding regions.

MACRO_HOW_TO_INSTALL

MACRO_PEDIGREE

MACRO_H2(Parameters)

MACRO_PARSE_CONFIG(../../../lib/nfconfigparser.jar, ../../../confs/by_workflow/wgselect.basic.cfg)

MACRO_EXECUTE(	--vcf /path/to/vcf \
	MACRO_COMMON_PARAMS \
	--pedigree "/path/to/file.ped" \
	--bed /path/to/input.bed 
)

MACRO_FOOTER
