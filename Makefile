all:
	$(MAKE) -C lib/nupack

clean:
	$(MAKE) -C lib/nupack clean
	rm -f nupyck/*.pyc
