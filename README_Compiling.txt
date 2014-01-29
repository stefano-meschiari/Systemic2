
Linux
=====

Installing from source
----------------------
(1), (3), (4) and (5) are available through your distribution's package manager (Sun's Java 6 or later, or OpenJDK 7 or later should be fine). Edit the Makefile to ensure all paths are set correctly (gcc, gfortran, libraries and includes -- particularly GSL). Run make reqs:

make -f Makefile.linux 

which will compile SWIFT and f2c in the local/ directory, download the needed R packages and build the systemic library (you should find libsystemic.so in the directory afterwards). If everything is set up correctly, you should be able to open the GUI by running

sh Systemic.sh

(you may modify this file to change java settings.)

Installing issues with GSL
--------------------------
You might not have permission to install GSL on the system, or get GSL-related compiling errors (e.g. error complaining about the "nsimplex2" symbol being  undefined). It might be advisable to build a local installation of the GSL library. A copy of the GSL 1.15 source package is included in this distribution. To build it, type

make -f Makefile.linux gsl

(the GSL library and includes will be installed in local/). Subsequently clean & make the main systemic library again:

make -f Makefile.linux clean
make -f Makefile.linux

The Makefile should automatically pick up the locally installed GSL library.
