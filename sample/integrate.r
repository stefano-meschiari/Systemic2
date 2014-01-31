# ======================================================================
# This script shows how to write a Systemic R script.
#
# Launch with
# /usr/bin/Rscript integrate.r
# ======================================================================

# Path to the systemic loader script, fix it to your install location
systemic.rpath <- "~/Projects/Systemic2/R/systemic.r"
# Load package
source(systemic.rpath, chdir=TRUE)
# This script will print and plot its result; set noisy to FALSE
# to suppress output.
noisy <- TRUE

# Create a new kernel object; this is the master object that
# contains data (stellar mass, planetary elements, RV data, etc.)
# and launches calculations.
k <- knew()

# Mstar = 1 Msun
k$mstar <- 1.
# 0 JD epoch
k$epoch <- 0.

# Add planets with orbital elements by hand...
k[] <- c(period = 4.617, mass = 0.689, ma = 314.71, ecc = 0.01, lop = 99.54)
k[] <- c(period = 241.33, mass = 1.9, ma = 36.98, ecc = 0.268, lop = 246.15)
k[] <- c(period = 1274.58, mass = 3.75, ma = 231.48, ecc = 0.259, lop = 252.92)

# ...or, say, read them from a text table. Every row is a planet, every column is an
# orbital element, in order period, mass, mean anomaly, eccentricity, long. of peri,
# inclination and node.
#
# kels(k) <- read.table("elementstable.txt")

nyears <- 1e5
if (noisy) { cat(sprintf("Starting integration for %e years (this might take a while)...\n", nyears)) }

# Start integration. int.method can be one of SWIFTRMVS, BULIRSCHSTOER or RK89; SWIFTRMVS is a
# fixed timestep integrator, and we need to specify a dt.
#
# The integration span can be specified with a single duration, in days...
i_1 <- kintegrate(k, times=nyears*365.25, int.method=SWIFTRMVS, dt=0.1*min(k[, 'period']))

# ... or as a vector of specific points at which to sample the integration
# i_1 <- kintegrate(k, times=c(100, 1000, 10000) * 365.25, int.method=SWIFTRMVS, dt=0.1*min(k[, 'period']))

# Print and plot results of integration? Set noisy to FALSE for silent execution
if (noisy) {
    print(i_1)

    # Saves plot as PDF
    pdf()
    plot(i_1)
    dev.off()
}

# i_1 is an object contains the results of the integration.
# i_1$els[[n]] (where n = 1..nplanets) are the orbital elements at each step
# of the integration.
# i_1$xyz[[n]] are the astrocentric cartesian coordinates at each step of the
# integration.
# i_1$rv[[n]] are the RVs that would be observed at each step.
#
# Save the results in nicely formatted text tables (one for each planet). The
# last column is the time at which the element snapshot was taken.

for (i in 1:k$nplanets) {
    fn <- sprintf("planet_%d.txt", i)
    write.fmatrix(i_1$els[[i]], file=fn)
}
