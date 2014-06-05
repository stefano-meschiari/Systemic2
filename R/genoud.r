default.ranges <- list()
default.ranges[[PER]] = c(0.1, 1e4)
default.ranges[[MASS]] = c(1e-4, 1e3)
default.ranges[[MA]] = c(-360, 360)
default.ranges[[ECC]] = c(0, 0.95)
default.ranges[[LOP]] = c(-360, 360)
default.ranges[[INC]] = c(-360, 360)
default.ranges[[NODE]] = c(-360, 360)
default.ranges[['par']] = c(-1e3, 1e3)

kminimize.domain <- function(k, log.period=FALSE, log.mass=FALSE) {
    indices <- kminimized.indices(k)
    
    return(t(sapply(1:k$nrpars, function(i) {
        r <- krange(k, indices[1, i], indices[2, i])

        if (log.period && indices[1, i] != -1 && indices[2, i] == PER)
            r <- log10(r)
        if (log.mass && indices[1, i] != -1 && indices[2, i] == MASS)
            r <- log10(r)
        
        
        if (is.nan(r[1]))
            if (indices[1, i] == -1) r[1] <- default.ranges[['par']][1] else r[1] <- default.ranges[[indices[2,i]]][1]
        if (is.nan(r[2]))
            if (indices[1, i] == -1) r[2] <- default.ranges[['par']][2] else r[2] <- default.ranges[[indices[2,i]]][2]
        
        return(r)
    })))
}




kminimize.genoud <- function(k, minimize.function='default', log.period=TRUE, log.mass=TRUE, max.generations=100, ...) {
    .require.library('rgenoud')
    .check_kernel(k)
    stopifnot(k$nplanets > 0)
    indices <- kminimized.indices(k)
    per.indices <- indices[1,] != -1 & indices[2,] == PER
    mass.indices <- indices[1,] != -1 & indices[2,] == MASS
    on.exit(kupdate(k, calculate=TRUE))
    
    kminimize.default <- function(v, ...) {
        if (log.period)
            v[per.indices] = 10^v[per.indices]
        if (log.mass)
            v[mass.indices] = 10^v[mass.indices]

        K_setMinimizedValues(k$h, v)
        K_calculate(k$h)

        return(K_getMinValue(k$h))
    }

    if (minimize.function == 'default')
        minimize.function = kminimize.default

    Domain <- kminimize.domain(k, log.period=log.period, log.mass=log.mass)
        
    v <- genoud(minimize.function, k$nrpars, Domains=Domain,  max.generations=max.generations, ...)
    
    return(invisible(v))
}

kminimize.de <- function(k, minimize.function='default', log.period=TRUE, log.mass=TRUE, population=50*k$nrpars,
                         max.iterations=1000, F = 0.5, CR = 1,
                         wait=10, check.function=NULL, ...) {
    .check_kernel(k)
    stopifnot(k$nplanets > 0)
    indices <- kminimized.indices(k)
    per.indices <- indices[1,] != -1 & indices[2,] == PER
    mass.indices <- indices[1,] != -1 & indices[2,] == MASS
    on.exit(kupdate(k, calculate=TRUE))
    
    Domain <- kminimize.domain(k, log.period=log.period, log.mass=log.mass)
    
    kminimize.default <- function(v, ...) {
        if (log.period)
            v[per.indices] = 10^v[per.indices]
        if (log.mass)
            v[mass.indices] = 10^v[mass.indices]

        K_setMinimizedValues(k$h, v)
        K_calculate(k$h)
        if (!is.null(check.function)) 
            if (!check.function(k))
                return(NaN);

        return(K_getMinValue(k$h))
    }

   

    if (minimize.function=='default')
        minimize.function <- kminimize.default
    
    x <- lapply(1:population, function(...) {
        return(runif(k$nrpars, Domain[,1], Domain[,2]))
    })


    
    f <- sapply(x, minimize.function)
    print(summary(f))
    pm <- sapply(x, function(p) return(p[c(1, 6)]))
    xrange <- range(pm[1,])
    yrange <- range(pm[2,])
    
    for (reps in 1:max.iterations) {
        x <- lapply(1:length(x), function(me) {
            repeat {
                idx <- sample(1:length(x), 3)
                if (! any(idx == me)) break
            }

            R <- sample(1:k$nrpars, 1)
            v <- sapply(1:k$nrpars, function(i) {
                if (runif(1) < CR || i == R)
                    return(x[[idx[1]]][i] + F * (x[[idx[2]]][i] - x[[idx[3]]][i]))
                else
                    return(x[[me]][i])
            })

            fnew <- minimize.function(v)
            if (fnew < f[me]) {
                f[me] <<- fnew
                return(v)
            }
            else
                return(x[[me]])
        })

        print(summary(f))
        pm <- sapply(x, function(p) return(p[c(1,6)]))
        plot(pm[1,], pm[2,])
        
    }    
}
