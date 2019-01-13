all: lib/nupack/Makefile
	$(MAKE) -C lib/nupack
	cp lib/nupack/nupack.so nupyck/
	cp -r lib/nupack/parameters nupyck/

lib/nupack/Makefile:
	git submodule update --init

clean:
	$(MAKE) -C lib/nupack clean
	rm -f nupyck/nupack.so
	rm -rf nupyck/parameters
