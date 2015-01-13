all: gaip bin

gaip:
	$(MAKE) --directory=gaip

bin:
	$(MAKE) --directory=src install

clean:
	$(MAKE) --directory=src clean
	$(MAKE) --directory=gaip clean
	-rm -f ./bin/*
	-find . -name "*.pyc" -delete
