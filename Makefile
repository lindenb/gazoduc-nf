BASH=/bin/bash
.PHONY:test tests clean


tests: test
test:
	cd tests && $(MAKE)

clean:
	rm -rf tests-output
