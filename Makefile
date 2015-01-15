all: gaip bin

gaip:
	$(MAKE) --directory=gaip

bin:
	$(MAKE) --directory=src install

docs:
	$(MAKE) --directory=docs html

clean:
	$(MAKE) --directory=src clean
	$(MAKE) --directory=gaip clean
	$(MAKE) --directory=docs clean
	-rm -f ./bin/*
	-find . -name "*.pyc" -delete

.PHONY: all clean gaip bin docs
