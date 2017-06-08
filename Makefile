all: lib/nupack/Makefile
	$(MAKE) -C lib/nupack

lib/nupack/Makefile:
	git submodule update --init

install:
	pip install -e .

uninstall:
	pip uninstall -y nupyck

clean:
	$(MAKE) -C lib/nupack clean
	rm -f nupyck/*.pyc
	rm -rf nupyck.egg-info
	pip uninstall -y nupyck
