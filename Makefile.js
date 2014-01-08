# Systemic 2 Makefile for Javascript
# 2014, Stefano Meschiari (http://www.stefanom.org)
#
# Compiling for Javascript requires at least emscripten 1.7.0 and llvm-3.2. Current "portable SDK" version is 1.5.6; activate "incoming" SDK to install 1.7.0.

# Steps:
# (1) Configure + compile gsl:
# emconfigure ./configure
# emmake make
# 

# If you have installed GSL & company in non-standard places, add them here:
LIBS = -L./deps/lib -L./private/javascript -L/opt/local/include
INCLUDES = -I./deps/include -I./private/javascript -I/usr/local/include -I/opt/local/include 
LIBNAMES = -lm -lc -lf2c 
INSTALL_DIR = /usr/local/include
PACKAGE_DIR = ~/Downloads
FFLAGS='-c -g -fPIC'


SYSTEMICLIVE_DIR = ~/Projects/SystemicLive/
JS_DIR=$(SYSTEMICLIVE_DIR)/js/

-include $(SYSTEMICLIVE_DIR)/Makefile_include
OPTIMIZED_FLAGS = -s LINKABLE=1 -O3 --closure 0 $(INCLUDES) $(LIBS) -DJAVASCRIPT -std=c99 -DJAVASCRIPT 
GSL_OBJECTS = private/javascript/f2c/*.o private/javascript/interpolation/*.o private/javascript/statistics/*.o private/javascript/blas/*.o private/javascript/cblas/*.o private/javascript/block/*.o private/javascript/err/*.o private/javascript/histogram/*.o private/javascript/matrix/*.o private/javascript/multifit/*.o private/javascript/multimin/*.o private/javascript/ode-initval2/*.o private/javascript/randist/*.o private/javascript/rng/*.o private/javascript/roots/*.o private/javascript/sort/*.o private/javascript/sys/*.o private/javascript/utils/*.o private/javascript/vector/*.o 

SYSFLAGS=$(OPTIMIZED_FLAGS)

#UPDATE = --update --java
UPDATE =

ALLOBJECTS = objects/periodogram.o objects/extras.o objects/mercury.o objects/integration.o objects/mcmc.o objects/utils.o objects/simplex.o objects/kernel.o objects/bootstrap.o objects/kl.o objects/qsortimp.o objects/lm.o objects/lm.o objects/ode.o objects/odex.o objects/sa.o objects/de.o

JS_FILES = ui help systemic

js: src/*.c src/*.h  $(ALLOBJECTS) support
	$(CC) $(CCFLAGS) $(SYSFLAGS) $(GSL_OBJECTS) $(EXP_FUNCTIONS) src/javascript.c  objects/*.o -o $(SYSTEMICLIVE_DIR)/js/libsystemic.html $(LIBS) $(LIBNAMES)  --embed-file datafiles 
	cd $(JS_DIR); sh minify.sh
	lua lua/minify_safely.lua $(JS_DIR)/libsystemic.js >$(JS_DIR)/libsystemic.min.js

support: 
	cd lua; luajit-2.0.0-beta9 jsparse.lua; luajit-2.0.0-beta9 list_sys.lua

test: src/*.c src/*.h $(ALLOBJECTS)
	$(CC) $(CCFLAGS) $(SYSFLAGS) $(GSL_OBJECTS) src/test.c -o packaging/forweb/systemic_test.js objects/*.o $(LIBS) $(LIBNAMES)  --embed-file datafiles --embed-file ./private/hd141399.fit

objects/kernel.o: src/kernel.c
	$(CC) $(CCFLAGS) $(SYSFLAGS) -c -o objects/kernel.o src/kernel.c 

objects/qsortimp.o: src/qsortimp.c
	$(CC) $(CCFLAGS) $(SYSFLAGS) -c -o objects/qsortimp.o src/qsortimp.c 

objects/extras.o: src/extras.c
	$(CC) $(CCFLAGS) $(SYSFLAGS) -c -o objects/extras.o src/extras.c 

objects/periodogram.o: src/periodogram.c
	$(CC) $(CCFLAGS) $(SYSFLAGS) -c -o objects/periodogram.o src/periodogram.c 

objects/simplex.o: src/simplex.c
	$(CC) $(CCFLAGS) $(SYSFLAGS) -c -o objects/simplex.o src/simplex.c 

objects/integration.o: src/integration.c
	$(CC) $(CCFLAGS) $(SYSFLAGS) -c -o objects/integration.o src/integration.c 

objects/mercury.o: src/mercury.c
	$(CC) $(CCFLAGS) $(SYSFLAGS) -c -o objects/mercury.o src/mercury.c

objects/utils.o: src/utils.c
	$(CC) $(CCFLAGS) $(SYSFLAGS) -c -o objects/utils.o src/utils.c

objects/bootstrap.o: src/bootstrap.c
	$(CC) $(CCFLAGS) $(SYSFLAGS) -c -o objects/bootstrap.o src/bootstrap.c

objects/kl.o: src/kl.c
	$(CC) $(CCFLAGS) $(SYSFLAGS) -c -o objects/kl.o src/kl.c

objects/mcmc.o: src/mcmc.c
	$(CC) $(CCFLAGS) $(SYSFLAGS) -c -o objects/mcmc.o src/mcmc.c

objects/lm.o: src/lm.c
	$(CC) $(CCFLAGS) $(SYSFLAGS) -c -o objects/lm.o src/lm.c

objects/ode.o: src/ode.c
	$(CC) $(CCFLAGS) $(SYSFLAGS) -c -o objects/ode.o src/ode.c

objects/odex.o: src/odex.c
	$(CC) $(CCFLAGS) $(SYSFLAGS) -c -o objects/odex.o src/odex.c

objects/sa.o: src/sa.c
	$(CC) $(CCFLAGS) $(SYSFLAGS) -c -o objects/sa.o src/sa.c

objects/de.o: src/de.c
	$(CC) $(CCFLAGS) $(SYSFLAGS) -c -o objects/de.o src/de.c

.PHONY: clean cleanreqs

#f2c: ./local/lib/libf2c.a
#	-echo "Making f2c"
#	cd local/f2c; emmake make -f makefile.u ; cp libf2c.a ../lib; cp f2c.h ../include

clean:
	rm -rf test objects/*.o 
