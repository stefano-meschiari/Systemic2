Dependencies
============
(1) Recent GNU Compiler Collection (GCC), 4.1 or better, including gfortran
(2) F2C (available here: http://www.netlib.org/f2c/)
(3) GNU Scientific Library (GSL)
(4) R 2.15 or better
(5) Sun Java 1.6 or better

Mac OS X
========

All-in-one binary
-----------------
An all-in-one binary (with all the requisite libraries installed in the folder) is available inside the mac/ folder. Double-click SystemicGui.jar to start it. If you'd rather build from sources, you can...

Installing from source
----------------------
You can download (1) from http://hpc.sourceforge.net already compiled. (3) can be installed through MacPorts or Fink, or from source. (2), (4) and (5) are bundled with the package. (6) is part of Mac OS X since Snow Leopard (Mac OS X 10.6).

Check that the Makefile has all the paths set up correctly (in particular, the paths to the libraries and includes so that GSL can be located correctly). To compile f2c, LuaJIT and Swift, type

make -f Makefile.osx reqs

which should build the requirements in the local/ subdirectory. Afterwards, type 

make -f Makefile.osx 

to compile the systemic library (you should see libsystemic.dylib now). If everything is set up correctly, you should be able to open the GUI by double-clicking on SystemicGUI.jar.

Linux
=====

Installing from source
----------------------
(1), (3), (4) and (5) are available through your distribution's package manager (Sun's Java 6 or later, or OpenJDK 7 or later should be fine). Edit the file Makefile.linux to ensure all paths are set correctly (gcc, gfortran, libraries and includes -- particularly GSL). Run make reqs:

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
