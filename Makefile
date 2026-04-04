BASH=/bin/bash
.PHONY:test tests clean doc schema


tests: test
test:
	cd tests && $(MAKE)

schema: src/xsl/metadata2json.xslt
	find ./workflows -type f -name "nextflow_schema.xml" |\
		while read F ; do xsltproc --xinclude src/xsl/metadata2json.xslt "$${F}"  | python3 -m json.tool > "$${F%.xml}.json" ; done

doc:
	source ${HOME}/.bash_profile && \
	micromamba activate NFCORE && \
	find ./workflows -type f -name "nextflow_schema.json"  -exec dirname '{}' ';' |\
		while read DIR ; do nf-core pipelines schema lint "$${DIR}/nextflow_schema.json" && nf-core pipelines schema docs --format markdown -o "$${DIR}/README.md" --force -c 'parameter,description,type,default,required,pattern' "$${DIR}" || true ; done


clean:
	rm -rf tests-output tests/.nextflow
