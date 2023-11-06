ifeq (${HOSTNAME},login-01.compute.bird2.prive)
JAVAC=/sandbox/apps/systeme/binaries/java/1.8.0_131/bin/javac
JAVACC=/LAB-DATA/BiRD/users/lindenbaum-p/packages/javacc/javacc.jar
else
JAVAC=javac
JAVACC=${HOME}/packages/javacc/target/javacc.jar
endif


../lib/nfconfigparser.jar : $(addsuffix .java,$(addprefix ./main/java/gazoduc/,config/NFConfigParser))
	mkdir -p $(dir $@)
	rm -rf TMP
	mkdir TMP
	$(JAVAC) -d TMP -sourcepath main/java $^
	jar cvf $@ -C TMP .
	rm -rf TMP

main/java/gazoduc/config/NFConfigParser.java : ./main/jj/ConfigParser.jj $(JAVACC)
	mkdir -p $(dir $@)
	java  -cp $(JAVACC) javacc -OUTPUT_DIRECTORY:$(dir $@) -JDK_VERSION:1.8 $<

