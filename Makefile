BASH=/bin/bash
.PHONY:test tests


tests:test
test:
	cd tests && $(MAKE) stub hard_tests
