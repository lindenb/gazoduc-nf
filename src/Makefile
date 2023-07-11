ifeq (${HOSTNAME},login-01.compute.bird2.prive)
JAVAC=/sandbox/apps/systeme/binaries/java/1.8.0_131/bin/javac
else
JAVAC=javac
endif

../doc/javadoc/index.html : ../lib/gazoduc.jar ./main/java/gazoduc/package-info.java
	mkdir -p $(dir $@)
	javadoc -notimestamp -d $(dir $@) -sourcepath main/java gazoduc

../lib/gazoduc.jar : $(addsuffix .java,$(addprefix ./main/java/gazoduc/,Gazoduc SAMSequenceDictionary SAMSequenceRecord config/NFConfigParser)) 
	mkdir -p $(dir $@)
	rm -rf TMP
	mkdir TMP
	$(JAVAC) -d TMP -sourcepath main/java $< main/java/gazoduc/config/NFConfigParser.java
	jar cvf $@ -C TMP .
	rm -rf TMP

main/java/gazoduc/config/NFConfigParser.java : ./main/jj/ConfigParser.jj
	mkdir -p $(dir $@)
	java  -cp /LAB-DATA/BiRD/users/lindenbaum-p/packages/javacc/javacc.jar javacc -OUTPUT_DIRECTORY:$(dir $@) -JDK_VERSION:1.8 $<
	
