source("~/Projects/Systemic2/R/systemic.r", chdir=TRUE)

# Set up a new kernel object
k <- knew()
k$mstar <- 0.95
k$epoch <- 0.

# Add planets
k[] <- c(period=38., mass=0.053, ma=214.07, ecc=0.22, lop=303.75, inc=90, node=0)
k[] <- c(period=122., mass=0.07, ma=33.57, ecc=0.36, lop=331.22, inc=90, node=0)

# Integrate for 10^5 years
nyears <- 1e4
#int <- kintegrate(k, nyears*365.25, int.method=SWIFTRMVS, dt=min(k[, 'period'])*0.05)
# plot(int)

# Grid points for planet #2: eccentricity between 0
# and 0.95, long. of pericenter between 0 and 360,
# sampled with 10 points respectively.
npoints <- 10
ecc <- seq(0, 0.95, length.out=npoints)
lop <- seq(0, 360, length=npoints)
grid <- expand.grid(ecc, lop)


# Return a list of 'survival times' for each 
# realization of the system with the modified orbital
# elements of planet #2
times <- mapply(function(e, l) {
    k2 <- kclone(k)
    k2[2, 'ecc'] <- e
    k2[2, 'lop'] <- l

    i <- kintegrate(k2, nyears*365.25, int.method=SWIFTRMVS, dt=min(k[, 'period'])*0.05)

    cat(sprintf("e = %.2f, lop = %.2f, survival time = %.2f\n", e, l, i$survival.time / 365.25))
    return(i$survival.time)
}, grid[[1]], grid[[2]])

# Cut the survival times by eccentricity...
times.by.ecc <- split(times, grid[[1]])

# ...and take the median over the phases
times.by.ecc.median <- sapply(times.by.ecc, median)

# Make a plot
plot(ecc, times.by.ecc.median / 365.25, type='l',
     xlab='Eccentricity', ylab='Survival time [yrs]', log='y')
