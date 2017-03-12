all: docs

gen-doc:
	docs/build-sphinx.sh

docs: gen-doc
	@$(MAKE) --directory=docs html

.PHONY: all docs
