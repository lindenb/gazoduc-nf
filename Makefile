BASH=/bin/bash
.PHONY:doc
define M4DOC
cd "$1" && m4 -P README.m4 > README.md && git add README.md
endef

all: doc

doc:
	$(call M4DOC,workflows/wgselect/basic)
