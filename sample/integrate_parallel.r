source("~/Projects/Systemic2/R/systemic.r", chdir=TRUE)
require(parallel)
# Set up a new kernel object
k <- knew()
k$mstar <- 0.95
k$epoch <- 0.

# Add planets
k[] <- c(period=38., mass=0.053, ma=214.07, ecc=0.219, lop=303.74, inc=90, node=0)
k[] <- c(period=122., mass=0.07, ma=33.57, ecc=0.36, lop=331.22, inc=90, node=0)

# Integrate for 5e4 years
nyears <- 5e4
int <- kintegrate(k, nyears*365.25, int.method=SWIFTRMVS, dt=min(k[, 'period'])*0.1)
pdf('integration_1e5.pdf')
plot(int)
dev.off()

# Grid points for planet #2: 
per <- seq(38, 122, length.out=40)
ma <- seq(0, 360, length=50)
grid <- expand.grid(per, ma)

# Return a list of 'survival times' for each 
# realization of the system with the modified orbital
# elements of planet #2
times <- mcmapply(function(p, ma) {
    k2 <- kclone(k)
    k2[2, 'period'] <- p
    k2[2, 'ma'] <- ma

    i <- kintegrate(k2, nyears*365.25, int.method=SWIFTRMVS, dt=min(k[, 'period'])*0.1)

    cat(sprintf("per = %.2f, ma = %.2f, survival time = %.2f\n", p, ma, i$survival.time / 365.25))
    return(i$survival.time)
}, grid[[1]], grid[[2]], mc.cores=8)

# Cut the survival times by period...
times.by.per <- split(times, grid[[1]])

# ...and take the median over the phases
times.by.per.median <- sapply(times.by.per, median)

# Save the results to file
save.image('results')

# Make a plot, using the default Systemic plot style
par(systemic.par)
plot(per, times.by.per.median / 365.25, type='l',
     xlab='Period [d]', ylab='Survival time [yrs]', log='y')

