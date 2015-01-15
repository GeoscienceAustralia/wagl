all: gaip bin

gaip:
	$(MAKE) --directory=gaip

bin:
	$(MAKE) --directory=src install

gen-doc:
	docs/build-sphinx.sh

docs: gen-doc
	$(MAKE) --directory=docs html

clean:
	$(MAKE) --directory=src clean
	$(MAKE) --directory=gaip clean
	$(MAKE) --directory=docs clean
	-rm -f ./bin/*
	-find . -name "*.pyc" -delete

.PHONY: all clean gaip bin docs
