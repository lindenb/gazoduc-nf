ifeq (${HOSTNAME},login-01.compute.bird2.prive)
JAVAC=/sandbox/apps/systeme/binaries/java/1.8.0_131/bin/javac
else
JAVAC=javac
endif

../doc/javadoc/index.html : ../lib/gazoduc.jar ./main/java/gazoduc/package-info.java
	mkdir -p $(dir $@)
	javadoc -notimestamp -d $(dir $@) -sourcepath main/java gazoduc

../lib/gazoduc.jar : $(addsuffix .java,$(addprefix ./main/java/gazoduc/,Gazoduc Genome SAMSequenceDictionary Genomes SAMSequenceRecord))
	mkdir -p $(dir $@)
	rm -rf TMP
	mkdir TMP
	$(JAVAC) -d TMP -sourcepath main/java $<
	jar cvf $@ -C TMP .
	rm -rf TMP
