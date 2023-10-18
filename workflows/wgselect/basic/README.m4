m4_define(`MACRO_WORKFLOW',wgselect)m4_dnl
m4_define(`MACRO_MAIN_NF',gazoduc-nf/workflows/wgselect/basic/wgselect.basic.nf)m4_dnl
m4_include(../../../src/main/m4/common.m4)m4_dnl
md_h1(MACRO_WORKFLOW)

md_h2(About)


HOW_TO_INSTALL

md_make

$ module purge
$ module load nextflow \
$ export NXF_ANSI_LOG=false
$ nextflow -c "../../gazoduc-nf/confs/login-01.compute.bird2.prive.cfg" -c ../../gazoduc-nf/confs/by_workflow/wgselect.basic.cfg run \
	-resume \
	-work-dir /SCRATCH-BIRD/users/lindenbaum-p/work/NEXTFLOW/2023/20230906.wgselect.hs37d5/work \
	/PATH/TO/MACRO_MAIN_NF \
	--vcf /LAB-DATA/BiRD/shares/ITX/u1087/lindenb/20211117.brs.mitral.chopin.genmed.gatk4.hs37d5/20211117.brs.mitral.chopin.genmed.gatk4.hs37d5.list \
	--pedigree "/SCRATCH-BIRD/users/lindenbaum-p/work/NEXTFLOW/2023/20230906.wgselect.hs37d5/work/20230906.wgselect.hs37d5.in.ped" \
	--genomeId hs37d5 \
	--prefix "20230906.wgselect.hs37d5." \
	--bed /SCRATCH-BIRD/users/lindenbaum-p/work/NEXTFLOW/2023/20230906.wgselect.hs37d5/work/scn5a.bed \
	--publishDir "/SCRATCH-BIRD/users/lindenbaum-p/work/NEXTFLOW/2023/20230906.wgselect.hs37d5/work"

