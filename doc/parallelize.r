# Example parallelized script.
#
# First load a kernel in the GUI
# with at least two planets (stop otherwise).
stopifnot(k$nplanets >= 2)

# Integrate for 10000 years
time <- 10000 * 365.25

# Parallelize over 8 threads
ncores <- 8

# Make a stability map over the first planet's period
# and mass. This is made by creating a list of kernels K
# with the respective periods and masses.

K <- list(.)

