m4_define(`md_hash', m4_changecom()#m4_changecom(#))m4_dnl
m4_define(`md_pre', m4_changequote([,])[m4_changequote([,])```[$1]m4_changequote(`,')]m4_changequote(`,'))m4_dnl
m4_define(`md_code', m4_changequote([,])[m4_changequote([,])`[$1]`m4_changequote(`,')]m4_changequote(`,'))m4_dnl
m4_define(`md_h1',m4_changecom()#m4_changecom(#) $1
)
m4_define(`md_h2',m4_changecom()##m4_changecom(#) $1
)
m4_define(`md_h3',m4_changecom()###m4_changecom(#) $1
)m4_dnl
m4_define(`md_quote',
> $1
)m4_dnl
m4_define(`md_france',:flag_fr:)m4_dnl
m4_define(`md_make',`md_pre(make)
$1
md_pre')m4_dnl
m4_define(`MACRO_HOW_TO_INSTALL',`
md_h2(Installation)

If it was not already done, clone the repo.
m4_pre(bash)
git clone https://gitlab.univ-nantes.fr/pierre.lindenbaum/gazoduc-nf.git
m4_pre

or, if you're using github

m4_pre(bash)
git clone https://github.com/lindenb/gazoduc-nf.git
m4_pre



')m4_dnl
