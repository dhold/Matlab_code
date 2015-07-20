include make.config


SOURCENAME=phi-source
BINARYNAME=phi-linux-x86_64

all:
	cd src;\
	make

doc:
	cd doc;\
	make

static:
	cd src;\
	make static

debug:
	cd src;\
	make debug

debug_install:
	cd src;\
	make debug;\
	echo "Copying phi to bin/";\
	cp phi ../bin	

install: 
	mkdir -p bin;\
	cd src;\
	make;\
	echo "Copying phi to bin/";\
	cp phi ../bin	

distrib: clean doc
	mkdir -p $(SOURCENAME)
	cp -r src $(SOURCENAME)
	cp -r doc $(SOURCENAME)
	cp -r examples $(SOURCENAME)
	cp license.txt README.txt notes.txt  $(SOURCENAME)
	cp Makefile make.config $(SOURCENAME)
	tar czf $(SOURCENAME).tar.gz $(SOURCENAME) --wildcards\
						       --exclude CVS\
						       --exclude $(BINARYNAME)\
						       --exclude $(SOURCENAME).tar.gz\
						       --exclude $(BINARYNAME).tar.gz
	rm -rf $(SOURCENAME)

binary_distrib: clean static doc
	mkdir -p $(BINARYNAME)
	mkdir -p $(BINARYNAME)/doc
	mkdir -p $(BINARYNAME)/bin
	cp -r src/phi $(BINARYNAME)/bin
	cp doc/phi_ug.pdf $(BINARYNAME)/doc
	cp -r examples $(BINARYNAME)
	cp license.txt README.txt notes.txt  $(BINARYNAME)
	cp Makefile make.config $(BINARYNAME)
	tar czf $(BINARYNAME).tar.gz $(BINARYNAME) --wildcards\
						       --exclude CVS\
						       --exclude $(SOURCENAME)\
						       --exclude $(SOURCENAME).tar.gz\
						       --exclude $(BINARYNAME).tar.gz
	rm -rf $(BINARYNAME)
	

clean:
	rm -f bin/phi
	rm -rf $(SOURCENAME)
	rm -rf $(BINARYNAME)
	cd src;\
	make clean

