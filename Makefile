all:
	$(MAKE) -C lib/nupack
	pip install -e .

clean:
	$(MAKE) -C lib/nupack clean
	rm -f nupyck/*.pyc
	pip uninstall -y nupyck
	rm -rf nupyck.egg-info
