# This Makefile is part of the NUPACK software suite
# Copyright (c) 2007 Caltech. All rights reserved.
# Coded by: Justin Bois 1/2007
#
# Makefile for making all NUPACK 3.0.5 subdirectories:
# src/thermo,  src/shared,  src/design


all: 	
	$(MAKE) -C src/
	gcc -fPIC -shared -g -o nupack.so \
  src/shared/*.o src/shared/rng/*.o src/thermo/utils/*.o

clean: 
	$(MAKE) -C lib/ clean
	$(MAKE) -C src/ clean
	rm -f nupack.so

