# in each subdirectory, look for a Makefile.  If found, run it.

DIRECTORIES = $(sort $(dir $(wildcard */)))

all: build
	for dir in $(DIRECTORIES) ; do \
		if [ -f ./$$dir/Makefile ]; then (echo "found Makefile in $$dir...\n"; cd $$dir; make all; echo "\n"); fi; \
	done

build:
	mkdir ./build

clean:
	for dir in $(DIRECTORIES) ; do \
		if [ -f ./$$dir/Makefile ]; then (echo "found Makefile in $$dir...\n"; cd $$dir; make clean; echo "\n"); fi; \
	done
	rm -rf ./build
