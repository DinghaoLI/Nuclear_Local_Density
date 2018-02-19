all: src doc

src:
	$(MAKE) -C src

doc:
	doxygen ./Doxyfile

test:
	$(MAKE) -C src test

.PHONY: clean src doc

clean:
	$(MAKE) -C src clean
	$(MAKE) -C doc clean

style:
	astyle --options=astyle.conf ./src/*.cpp ./src/*.h