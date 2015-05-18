## Creating/opening/saving kernel objects
* [knew](#knew) - Creates a new kernel object. 
* [kclone](#kclone) - Returns a new, independent copy of the kernel, with the same data and planets loaded. 
* [kload.old](#kload.old) - Loads an old-style fit (from the previous version of Systemic). 
* [kload.datafile](#kload.datafile) - Loads the datafiles contained in a .sys file (see, e.g., the .sys files contained in the datafiles folder). 
* [kload](#kload) - Loads a kernel (previously saved with [ksave](#ksave)) from disk. 
* [ksave](#ksave) - Saves a kernel (or a list of kernels) to a file. 

## Setting fit parameters (planets, offsets, etc.)
* [kadd.planet](#kadd.planet) - Adds a new body with the given elements (you can specify period, mass, ecc, lop, inc, node) 
* [kremove.planet](#kremove.planet) - Removes the idx-th planet 
* [kels](#kels) - Returns a matrix of orbital elements 
* [kels<-](#kels<-) - Sets the orbital elements of the kernel 
* [kernel\[idx1, idx2\]](#kernel\[idx1, idx2\]) - Subsetting a kernel object with brackets returns the values of the orbital elements. 
* [k\[idx1, idx2\] <- value](#k\[idx1, idx2\] <- value) - The parameters of a kernel object can be set by subsetting it with brackets and assigning values.  
* [kallels](#kallels) - Returns a matrix of orbital elements, including derived orbital elements 
* [kpars](#kpars) - Returns a vector of parameters 
* [kflag](#kflag) - Returns/sets the flag (active, minimized or inactive) for an orbital element or parameter. 
* [krange](#krange) - Returns the allowed range of the given parameter. 
* [krange<-](#krange<-) - Sets the allowed range of the given parameter. 
* [keltype](#keltype) - Sets the orbital elements format. 

## Minimization
* [kminimize](#kminimize) - Minimizes the chi^2 of the fit. 
* [kminimize1d](#kminimize1d) - Minimizes the specified parameter 

## Error estimation
* [kcrossval.l1o](#kcrossval.l1o) - Runs the "leave-one-out" cross validation algorithm. 
* [kbootstrap](#kbootstrap) - Runs the bootstrap routine on the given kernel. 
* [kmcmc](#kmcmc) - Runs the MCMC routine on the given kernel. 

## Loading and manipulating data
* [kadd.data](#kadd.data) - Adds a new dataset to the kernel. 
* [kdata](#kdata) - Returns a matrix containing the idx-th dataset (or all the data) 
* [kdata<-](#kdata<-) - Sets a matrix containing the idx-th dataset (or all the data) 

## Fit parameters (chi^2, jitter, etc.)
* [k\$property](#k\$property) - Use the $ operator to access the following properties of the kernel. 
* [k\$property <- value](#k\$property <- value) - Use the $ operator to set the following properties of the kernel. 
* [kcalculate](#kcalculate) - Recalculates or updates the statistics for the kernel (e.g. chi^2, rms, jitter, radial velocity response, etc.). 

## Periodogram
* [kperiodogram](#kperiodogram) - Returns a periodogram of the supplied time series. 
* [kperiodogram.boot](#kperiodogram.boot) - Returns a periodogram of the supplied time series, where the false alarm probabilities are estimated using a bootstrap method. 

## Constants

## Other
* [kstep](#kstep) - Returns the "step" for a given parameter.
* [kstep<-](#kstep<-) - Sets the "step" for a given parameter.
* [kselect](#kselect) - Selects (make available for minimization) the given parameters.
* [kdeselect](#kdeselect) - Deselects (exclude from minimization) the given parameters.
* [kxyz](#kxyz) - Returns the cartesian coordinates of the bodies in the system. 
* [krvcurve](#krvcurve) - Calculates the radial velocity curve over the specified time vector.

<hr>
<a name='knew'></a>

## knew
**knew() **

Creates a new kernel object.

### Returns:

A new kernel object. A kernel object contains both the fit parameters and the data.

<hr>
<a name='kadd.planet'></a>

## kadd.planet
**kadd.planet(k, period = 300, mass = 1, ma = 0, ecc = 0, lop = 0,
                       inc = 90, node = 0, tperi = NA, K = NA, a = NA) **

Adds a new body with the given elements (you can specify period, mass, ecc, lop, inc, node)

Alternative syntax:
k[] <- c(period = 300, mass = 1, ...)

### Arguments:


- k: the kernel to add the planet to
- period: period of planet (days)

- mass: mass of planet (Jupiter masses)
- ma: mean anomaly (deg)

- ecc: eccentricity
- lop: longitude of pericenter (deg)

- inc: inclination (deg)
- node: node (deg)

- tperi: time of passage through pericenter (days)
- K: semi-amplitude (m/s)

- a: semi-major axis (AU)

<hr>
<a name='kremove.planet'></a>

## kremove.planet
**kremove.planet(k, idx) **

Removes the idx-th planet
Alternative syntax: k[idx] <- NULL

### Arguments:


- k: kernel to remove the planet from
- idx: index of the planet to remove (starting at 1)

<hr>
<a name='kels'></a>

## kels
**kels(k, keep.first = FALSE) **

Returns a matrix of orbital elements
Alternative syntax: k[]

### Arguments:

k: kernel to read the orbital elements from

### Returns:

A matrix of orbital elements; each row is a planet's orbital elements.

<hr>
<a name='kels<-'></a>

## kels<-
**`kels`(k) <-  value**

Sets the orbital elements of the kernel

### Arguments:


- k: kernel where the orbital elements should be set
- elements: a matrix of orbital elements

<hr>
<a name='kernel\[idx1, idx2\]'></a>

## kernel\[idx1, idx2\]
**kernel[idx1, idx2]**

Subsetting a kernel object with brackets returns the values of the orbital elements.

* k[] returns a matrix of elements (like kels)
* k[n, ] returns the elements for the n-th planet

* k[, m] returns the m-th element for all planets
* k[n, m] returns the m-th element for the n-th planet

* k['par'] returns a vector of all the parameters (e.g. offsets, linear trend, etc.)
* k['par', j] returns the j-th parameter

<hr>
<a name='k\[idx1, idx2\] <- value'></a>

## k\[idx1, idx2\] <- value
**k[idx1, idx2] <- value**

The parameters of a kernel object can be set by subsetting it with brackets and assigning values.

- k[] <- matrix sets the orbital elements to the specified matrix
- k[] <- c(period=..., mass=..., ...) adds a new planet with the specified orbital elements

- k[n] <- NULL removes the n-th planet
- k[n, m] <- v sets the m-th element for the n-th planet to v

- k[, m] <- v sets the m-th element to v for all planets
- k['par'] <- v sets parameter values to the vector v

- k['par', n] <- v sets the n-th parameter to v

<hr>
<a name='k\$property'></a>

## k\$property
**k$property**

Use the \$ operator to access the following properties of the kernel.
Read-only properties:

* k\$nplanets	Number of planets
* k\$chi2 		Current reduced chi^2 value (normalized by (k\$ndata - k\$nrpars))

* k\$chi2nr		Non-reduced chi^2
* k\$rms		Current RMS value

* k\$jitter	Current jitter value
* k\$loglik	Current log likelihood (multiplied by -1)

* k\$ks.pvalue	Current p-value of the KS test comparing normalized residuals to a unit gaussian
* k\$ndata		Number of data points

* k\$nrvs		Number of RV data points
* k\$ntts		Number of central transits

* k\$nsets		Number of data sets loaded
* k\$chi2rvs	Chi^2 (RVs only)

* k\$chi2tts	Chi^2 (central transits only)
* k\$trange	Time range of the compiled dataset

* k\$epoch Epoch of the fit (JD)
* k\$mstar Mass of the star (Msun)

* k\$int.method	Integration method (possible values: KEPLER, RK45, RK89)
* k\$element.type Coordinate type (possible values: ASTROCENTRIC, JACOBI)

* k\$min.func Function minimized by [kminimize.](#kminimize.) Possible values are "chi2" (default), "rms", or a function that takes a kernel and returns a number.
* k\$nrpars	"Degrees of freedom" parameter used to calculate reduced chi^2. It is equal to the number of all the parameters that are marked as ACTIVE or MINIMIZE

<hr>
<a name='k\$property <- value'></a>

## k\$property <- value
**k$property <- value**

Use the \$ operator to set the following properties of the kernel.
Settable properties:

* k\$int.method	Integration method (possible values: KEPLER, RK45, RK89)
* k\$element.type Coordinate type (possible values: ASTROCENTRIC, JACOBI)

* k\$epoch		Epoch in JD
* k\$mstar		Mass of the star in solar masses

* k\$minimize.func	Function or property string (one of chi2, rms, jitter, chi2nr, chi2rvs, chi2tts, loglik) which specifies the value to minimize
* k\$min.func Function minimized by [kminimize.](#kminimize.) Possible values are "chi2" (default), "rms", or a function that takes a kernel and returns a number.

<hr>
<a name='kallels'></a>

## kallels
**kallels(k, keep.first = F) **

Returns a matrix of orbital elements, including derived orbital elements

### Arguments:


- k: kernel

### Returns:

A matrix with all orbital elements (including semi-major axis, K, etc.)

<hr>
<a name='kpars'></a>

## kpars
**kpars(k) **

Returns a vector of parameters

### Arguments:


- k: kernel

### Returns:

A vector of parameters (e.g. offsets, trends)

<hr>
<a name='kadd.data'></a>

## kadd.data
**kadd.data(k, data, type = NA) **

Adds a new dataset to the kernel.

### Arguments:


- k: kernel to add data to
- data: either the path to a data file in textual format, or a matrix containing data

- type: type of data contained in the data argument. One of the RV or TIMING constants

<hr>
<a name='kdata'></a>

## kdata
**kdata(k, idx = 'all') **

Returns a matrix containing the idx-th dataset (or all the data).
The columns are:

- TIME (time of observation)
- VAL (measurement)

- ERR (uncertainty on the measurement)
- SVAL (measurement, shifted vertically by the offset parameters)

- PRED (value of the fit at the time of measurement)
- SET (index of the dataset that contains the data point)

### Arguments:


- k: kernel
- idx: either the index of the dataset you are interested in, or 'all' to get all the data compiled from the loaded datasets.

<hr>
<a name='kdata<-'></a>

## kdata<-
**`kdata`(k, idx="all") <-  value**

Sets a matrix containing the idx-th dataset (or all the data). See also [kdata.](#kdata.)

### Arguments:


- k: kernel
- idx: either the index of the dataset you are interested in, or 'all' to get all the data compiled from the loaded datasets.

- value: a matrix of observations.

<hr>
<a name='kcalculate'></a>

## kcalculate
**kcalculate(k) **

Recalculates or updates the statistics for the kernel (e.g. chi^2, rms, jitter, radial velocity response, etc.).

### Arguments:


- k: kernel

<hr>
<a name='kclone'></a>

## kclone
**kclone(k) **

Returns a new, independent copy of the kernel, with the same data and planets loaded.

### Arguments:

k: a kernel

### Returns:

A new kernel object that clones the input kernel object.

<hr>
<a name='kload.old'></a>

## kload.old
**kload.old(file, datafiles.dir=paste(dirname(file), "/datafiles")) **

Loads an old-style fit (from the previous version of Systemic).

### Arguments:

file: path to the fit to be loaded
datafiles.dir: path to the datafiles (.sys and .vels files)

### Returns:

A new kernel object with data and parameters loaded according to the contents of file
InitialEpoch:")) {
InitialEpoch:")

<hr>
<a name='kload.datafile'></a>

## kload.datafile
**kload.datafile(k, datafile) **

Loads the datafiles contained in a .sys file (see, e.g., the .sys files contained in the datafiles folder).

### Arguments:


- k: kernel to load data into
- datafile: path to the .sys file

<hr>
<a name='kload'></a>

## kload
**kload(file, skip = 0, chdir=TRUE) **

Loads a kernel (previously saved with [ksave](#ksave)) from disk.

### Arguments:


- file: path to the file
- skip: if multiple kernels are saved in a single file, the index of the kernel to read (usually only one kernel is saved per file, so specify 0)

<hr>
<a name='ksave'></a>

## ksave
**ksave(k, file) **

Saves a kernel (or a list of kernels) to a file.

### Arguments:


- k: a single kernel, or a list of kernels to save
- file: file to save to. Multiple kernels might be saved to the same file.

<hr>
<a name='kperiodogram'></a>

## kperiodogram
**kperiodogram(k, per_type = "all", samples = getOption("systemic.psamples", 50000), pmin = getOption("systemic.pmin", 0.5), pmax = getOption("systemic.pmax", 1e4), data.flag = T_RV, timing.planet = NULL, val.col = SVAL, time.col = TIME, err.col = ERR, pred.col = PRED, plot = FALSE, print = FALSE, peaks = 25,
                         overplot.window=TRUE, .keep.h = FALSE) **

Returns a periodogram of the supplied time series.
If the first parameter is a kernel, then this function will return
periodogram of the loaded datasets (if per_type = "all") or the
residuals (if per_type = "res"). If the first parameter is a matrix,
then this function will return a periodogram of the columns of the
matrix.
The periodogram function will need three columns of data: a timestamp
column, a value column (e.g. the RV amplitude at that timestamp), and
an uncertainty column. The default indexes for those columns are TIME
[1], SVAL and ERR.
The false alarm probability returned is only an analytical
estimate. For more accurate false alarm probabilities, use
kperiodogram.boot (which uses a bootstrap method to estimate the
false alarm probabilities).

### Arguments:


- k: either a kernel object, or a matrix
- per_type: if k is a kernel object, one of "all" (periodogram of
the full compiled dataset) or "res" (periodogram of the
residuals)

- samples: number of periods (frequencies) at which to sample the
periodogram

- peaks: identifies the N tallest peaks in the periodogram
- pmin: minimum period at which to sample the periodogram

- pmax: maximum period at which to sample the periodogram
- data.flag: type of data (T_RV or T_TIMING)

- timing.planet: if the data is of type T_TIMING, specifies which
transits to use to calculate the periodogram

- val.col: the column to consider as the "value" column (by default, the SVAL column)
- time.col: the column to consider as the "time" column (by default, the TIME column)

- err.col: the column to consider as the "uncertainty" column (by default, the ERR column)
- pred.col: the column to consider as the "model value" column
(used if per_type = "res") (by default, the PRED column)

- plot: if TRUE, plot the periodogram after the calculation
- overplot.window: if TRUE, overplot the periodogram of the sampling
window

- print: if TRUE, pretty-prints the periodogram sorted by power.


### Returns:

A matrix with columns containing, respectively: period, power
at that period, (analytical) false alarm probability,
unnormalized power, tau, power of the sampling window at that
period

<hr>
<a name='kperiodogram.boot'></a>

## kperiodogram.boot
**kperiodogram.boot(k, per_type = "all", trials = 1e5, samples = getOption("systemic.psamples", 50000), pmin = getOption("systemic.pmin", 0.5), pmax = getOption("systemic.pmax", 1e4), data.flag = T_RV, timing.planet = NULL, val.col = SVAL, time.col = TIME, err.col = ERR, seed = 1, plot = FALSE, print = FALSE,
                              overplot.window=TRUE, peaks=25) **

Returns a periodogram of the supplied time series, where the false alarm probabilities are estimated using a bootstrap method.
If the first parameter is a kernel, then this function will return
periodogram of the loaded datasets (if per_type = "all") or the
residuals (if per_type = "res"). If the first parameter is a matrix,
then this function will return a periodogram of the columns of the
matrix.
The periodogram function will need three columns of data: a timestamp
column, a value column (e.g. the RV amplitude at that timestamp), and
an uncertainty column. The default indexes for those columns are TIME
[1], SVAL and ERR.
The false alarm probabilities are estimated by computing "trials"
periodograms of gaussian noise. The routine is automatically
parallelized.

### Arguments:

The arguments are the same as the [kperiodogram](#kperiodogram) function, plus the following:

- trials: number of periodograms to use in the boostrap estimation
- samples: number of periods (frequencies) at which to sample the
periodogram

### Returns:

A matrix with columns containing, respectively: period, power
at that period, false alarm probability (calculated in the bootstrap procedure),
unnormalized power, tau, power of the sampling window at that
period

<hr>
<a name='kflag'></a>

## kflag
**kflag(k, row, column) **

Returns/sets the flag (active, minimized or inactive) for an orbital element or parameter.
This corresponds to clicking one of the "semaphore" buttons next to the fit parameters in the user interface.
The value set can be one of INACTIVE (not minimized over, not counted
as a parameter in reduced chi^2), ACTIVE (not minimized over, counted
as a parameter in reduced chi^2), MINIMIZE (minimized over). By
default, many orbital elements are ACTIVE + MINIMIZE (these appear
as green in the GUI).

### Arguments:


- k: kernel object
- row: either the planet index, or "par" to specify a parameter
(e.g. an RV offset)

- col: either the orbital element index (one of `constants`)
- value: one of INACTIVE, ACTIVE, MINIMIZE, or ACTIVE + MINIMIZE

<hr>
<a name='kminimize'></a>

## kminimize
**kminimize(k, iters = 5000, algo = NA, de.CR = 0.2,
                      de.NPfac = 10, de.Fmin = 0.5, de.Fmax = 1.0, de.use.steps = FALSE,
                      sa.T0 = k$chi2, sa.alpha=2, sa.auto=TRUE, sa.chains=4) **

Minimizes the chi^2 of the fit.
kminimize uses one of the built-in algorithms to minimize the
reduced chi^2 of the fit. The algorithm is one of the following:

- SIMPLEX uses the Nelder-Meade algorithm (as implemented in GSL)
to search for a local minimum.

- LM uses the Levenberg-Marquardt algorithm (as implemented in GSL)
to search for a local minimum.

- SA uses a simple implementation of the simulated annealing
algorithm.

- DE uses a simple implementation of the differential evolution
algorithm.
The minimization algorithms may use the parameter steps set by
[kstep](#kstep) as initial scale parameters to explore the chi^2 landscape.
The target function to be minimized is defined by k\$min.func, which
has a value of "chi2" by default.

### Arguments:


- k: kernel to minimize
- iters: maximum number of iterations

- algo: one of SIMPLEX, LM, SA or DE. If none is specified, uses
the value in k\$min.method

- sa.T0: for SA, the initial temperature of the annealer
- sa.alpha: the index of the annealer (T = T0 (1 - (n/N)^alpha))

- sa.auto: automatically derive steps that produce a variation of chi^2 = 10% T0
- de.CR: crossover probability for DE

- de.Fmin, de.Fmax: differential weight for DE

<hr>
<a name='kcrossval.l1o'></a>

## kcrossval.l1o
**kcrossval.l1o(k, iters = 5000, algo = NA, type=NA) **

Runs the "leave-one-out" cross validation algorithm.

### Arguments:


- k: the kernel to run the routine on.
- iters: number of iterations (see [kminimize](#kminimize))

- algo: minimization algorithm (see @kminimize)
- type: if NA, returns the result of the cross-validation algorithm; otherwise, runs the cross-validation algorithm by removing planets of progressively smaller k.
", .k\$nplanets, " planets\n")
", .k\$nplanets, " planets\n")

<hr>
<a name='kminimize1d'></a>

## kminimize1d
**kminimize1d(k, row, column, algo = NA, iters = 5000) **

Minimizes the specified parameter

### Arguments:


- k: kernel
- row: either the index of the planet, or "par" to minimize a data parameter

- column: the planet parameter (e.g. 'period', 'ma', 'ecc') or the data parameter
- algo, iters: see [kminimize.](#kminimize.)

<hr>
<a name='krange'></a>

## krange
**krange(k, row, column) **

Returns the allowed range of the given parameter.
See [krange](#krange)<-.

<hr>
<a name='krange<-'></a>

## krange<-
**`krange`(k, row, column) <-  value**

Sets the allowed range of the given parameter.
A parameter range forces a parameter to lie within a specified
interval. The interval is specified as a two-element vector,
c(min, max). If either (or both) min or max is NaN, then
the minimum/maximum is not enforced.

### Arguments:


- k: kernel
- row: Either the index of the planet, or 'par' to specify the range of a data parameter.

- column: Either the name of the planet parameter (e.g. 'period', 'ma', etc.) or the name of the data parameter.
- value: A two-element vector c(min, max). min and max can be NaN.

<hr>
<a name='kstep'></a>

## kstep
**kstep(k, row, column) **

Returns the "step" for a given parameter.
See @`kstep<-`.

<hr>
<a name='kstep<-'></a>

## kstep<-
**`kstep`(k, row, column) <-  value**

Sets the "step" for a given parameter.
The "step" for a parameter specifies the typical step taken by
minimization routines for the parameter. For the SIMPLEX algorithm,
it represents the initial size of the simplex; for the LM algorithm,
it specifies the step over which gradients are calculated.

### Arguments:


- k: kernel
- row: either the planet index, or 'par' to specify the step of a data parameter.

- column: Either the name of the planet parameter (e.g. 'period', 'ma', etc.) or the name of the data parameter.
- value: The value of the step.

<hr>
<a name='kselect'></a>

## kselect
**kselect(k, row = "all", column = "all") **

Selects (make available for minimization) the given parameters.
This corresponds to the action of clicking on the semaphore buttons
in the user interface to select a parameter for minimization.

### Arguments:


- k: the kernel
- row: either the planet index, or 'par' to specify a data parameter.

- column: either the name of the planet parameter (e.g. 'period', 'ma', etc.) or the name of the data parameter.

<hr>
<a name='kdeselect'></a>

## kdeselect
**kdeselect(k, row = "all", column = "all") **

Deselects (exclude from minimization) the given parameters.
See [kselect.](#kselect.)

<hr>
<a name='kbootstrap'></a>

## kbootstrap
**kbootstrap(k, algo = NA, trials = 5000, warmup = 0, min_iter = 2000, plot = FALSE, print = FALSE, save=NA) **

Runs the bootstrap routine on the given kernel.
This function runs the bootstrap algorithm to estimate the
uncertainty on the parameters.

### Arguments:


- k: the kernel to run bootstrap on
- algo: the algorithm to use in minimization passes (see [kminimize](#kminimize))

- trials: the number of resampling trials
- plot: plots the resulting uncertainty object

- print: prints the resulting uncertainty object

<hr>
<a name='kmcmc'></a>

## kmcmc
**kmcmc(k, chains= 2, temps = 1, start = "perturb", noise=TRUE, skip.first = 1000, discard = k$nrpars * 10, R.stop = 1.1, 
                  min.length = 5000, max.iters = -1, auto.steps = TRUE, acc.ratio = 0.44, plot = FALSE, print = FALSE, save=NA,
                  debug.verbose.level = 1, random.log=TRUE) **

Runs the MCMC routine on the given kernel.
This function runs a simple implementation of MCMC on the kernel
to estimate uncertainty on the parameters.

### Arguments:


- k: the kernel to run bootstrap on, or a list of kernels with different parameters which represent the starting initial conditions.
- chains: number of chains to run in parallel

- skip.first: discard the first iterations
- R.stop: the Gelman-Rubin statistic used to estimate when to stop the routine.

- discard: only retain every n-th element of the chain
- min.length: minimum number of iterations

- acc.ratio: the acceptance ratio
- print: prints the resulting uncertainty object

- plot: plots the resulting uncertainty object

<hr>
<a name='kxyz'></a>

## kxyz
**kxyz(k) **

Returns the cartesian coordinates of the bodies in the system.
Coordinates are returned in internal units (Msun, AU, and day)

<hr>
<a name='keltype'></a>

## keltype
**keltype(k, value) **

Sets the orbital elements format.

### Arguments:


- k: kernel
- value: one of ASTROCENTRIC for astrocentric elements, or JACOBI for Jacobi elements.

<hr>
<a name='krvcurve'></a>

## krvcurve
**krvcurve(k, times=seq(from=min(ktrange(k)), to=max(ktrange(k)), length.out=getOption("systemic.rvsamples", 5000))) **

Calculates the radial velocity curve over the specified time vector.

### Arguments:


- k: the kernel
- times: a vector of times where to sample the radial velocity curve.

