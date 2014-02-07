Installation instructions for Systemic2
======================================= 
(c) 2012-2014, Stefano Meschiari http://www.stefanom.org

Note: if you have downloaded the Mac distribution of Systemic, you do
*not* need to compile the source code (it is already compiled for
you).


Introduction 
============ 
Systemic Console 2 (http://www.stefanom.org/systemic) is a new
software package for fitting and analyzing a variety of exoplanetary
data. It consists of a fast, parallelized C library (which you can
easily integrate into your own project), an R package providing a
high-level interface to the library, and an easy-to-use ”studio” user
interface for interactive use and plotting.

Please move the Systemic directory anywhere on your hard drive
(typically, somewhere inside your home folder) and follow the
instructions below to install Systemic.



Installing Systemic2 
==================== 
To install Systemic2, you will have to first compile its source
code. Please follow these steps:

(1) Please ensure that the following dependencies (i.e. pre-requisite
software that needs to be installed before Systemic) are met. You
should be able to find these packages in your distribution's package
manager.

*Dependencies*

- Recent GNU Compiler Collection (GCC), 4.1 or better, including
  gfortran. (Notice: clang will not work, for now) 
- F2C (available here: http://www.netlib.org/f2c/) 
- GNU Scientific Library (GSL), 1.15 or better (current version: 1.16) 
- R 2.15 or better (current version: 3.02) 
- Sun Java 1.6 or better (current version: 1.8)

Systemic will *not* work if these prerequisites are not met.

(2) Once you have verified that these dependencies are installed
correctly, edit the Makefile.linux file and modify these paths, if
necessary. Usually the default paths work, so you should be able to
skip this step.

- CC [the path of gcc] 
- FORTRAN [the path of gfortran] 
- R [the path of R] 
- LIBS [add the path where GSL and other system libraries are
  installed, if they are not already in the search path; usually
  /usr/local/lib or /opt/local/lib] 
- INCLUDES [add the path where GSL  and other include files are
  installed, if they are not already in the search path; usually 
  /usr/local/include or /opt/local/include]

(3) Cd to the Systemic directory and execute the Makefile, like so:

> make -f Makefile.linux

This will (a) unpack and compile some fortran dependencies, (b)
compile the Systemic source code, and (c) download and install some R
packages that are needed to use Systemic (an active Internet
connection is required).



Running Systemic 
================

Starting the user interface 
--------------------------- 
Start the Systemic point-and-click user interface by cd'ing into the
Systemic directory and typing:

> sh Systemic.sh

Running Systemic as an R package 
--------------------------------
Make sure to update the library search path to include the Systemic
directory (so that other programs, including R, can find 
libsystemic.so). You can do that by typing

> export LD_LIBRARY_PATH="/path/to/Systemic:$LD_LIBRARY_PATH"

(substitute the correct path to Systemic). You can add this line
to the shell startup file (e.g. ~/.bash_profile or ~/.bashrc).

Start R and type

> source('/path/to/Systemic/R/systemic.r', chdir=TRUE)

This will import all the Systemic routines (these routines usually
start with the 'k' letter).

Using Systemic as a C library 
----------------------------- 
You can use Systemic as a C library, so to use its algorithms from
within your project.  Headers are available in src (e.g.: systemic.h,
kernel.h, etc.). Most functions are replicated in the R package.

The Systemic library is called 'libsystemic.so' and is placed in the
Systemic directory.



Possible issues 
===============

Makefile complaining about paths 
-------------------------------- 
Make sure that you have set paths correctly in step (2).

SystemicGui.jar is crashing when double-clicked
----------------------------------------------- 
Use Systemic.sh to start Systemic. This script attempts to set up all
search paths correctly.

GSL 
--- 
You might not have permission to install GSL on the system, or
get GSL-related compiling errors (e.g. error complaining about the
"nsimplex2" symbol being undefined). It might be advisable to build a
local installation of the GSL library. A copy of the GSL 1.15 source
package is included in this distribution. To build it, type

make -f Makefile.linux gsl

(the GSL library and includes will be installed in
local/). Subsequently clean & make the main systemic library again:

make -f Makefile.linux clean make -f Makefile.linux

The Makefile should automatically pick up the locally installed GSL
library.

Systemic complains that it cannot find 'libsystemic.so'
-------------------------------------------------------
The Makefile might have failed at some point and did not produce
the compiled shared library needed by the user interface and R
package.

Try running the Makefile again, piping the output and scanning it 
for errors.


Need help?
==========
- File a bug at https://github.com/stefano-meschiari/Systemic2/issues
- Contact me; you can find my contact information at 
  http://www.stefanom.org



Enjoy! :)
